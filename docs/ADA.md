# Working on Ada (IIIT HPC)

What we verified on Ada on 12 Oct 2026, plus the relevant parts of the official
[Ada User Guide](https://hpc.iiit.ac.in/wiki/index.php?title=Ada_User_Guide). **The guide is partly out of date**
(old hostname, old partitions, Ubuntu 18.04 modules): trust this file where they differ. Needs the IIIT VPN off campus.

## Verified facts (smoke test `jobs/check_env.sbatch`, job 17358)

| Item | Value |
|---|---|
| SSH gateways | `ada-gw1.iiit.ac.in` (10.4.24.30), `ada-gw2.iiit.ac.in` (10.4.24.29); gateways run Rocky Linux 9.8. `ada.iiit.ac.in` refuses SSH |
| Home | `/home2/$USER` (NFS, shared by both gateways and compute nodes); `df` shows a 50 GB quota |
| GPU partition | **`u22`** (default): 77 nodes, Ubuntu 22.04.5. `short` / `long` from the guide no longer exist |
| RTX 2080 Ti nodes | `--constraint=2080ti` (45 nodes, mostly 4 GPUs each); other u22 nodes have 2–4 GPUs without a feature tag |
| GPU seen in job | NVIDIA GeForce RTX 2080 Ti, 11264 MiB, driver 570.211 (supports CUDA 12.x wheels) |
| Other partitions | `u22-cpu` (6 CPU nodes); `ihub` (3080 Ti), `plafnet2`/`rrc` (3090) are group-restricted |
| Modules (u22) | `u22/cuda/11.7, 11.8, 12.1, 12.4, 12.6, 12.9`; `u22/cudnn/9.1.0-cuda-12`, `9.7.0-cuda-12.8`; `u22/python/3.12.4`, `3.12.15`; system `python3` = 3.10.12 |
| Internet from compute nodes | **yes**: GDC, Hugging Face and PyPI all return HTTP 200 |
| Node-local disks | `/scratch` 1.8 TB, `/ssd_scratch` 880 GB (shared with other users' jobs on that node) |
| `/share1` | mounted on the gateway, **not** on compute nodes |
| Account / QoS | `research` / `medium`: **40 CPUs, 4 GPUs, 125 GB RAM per user across all jobs**, 20 jobs running/submitted |

## Connecting

`~/.ssh/config` on the laptop:

```
Host ada ada-gw
    HostName ada-gw1.iiit.ac.in
    User <ada-username>
    IdentityFile ~/.ssh/ada_claude
    IdentitiesOnly yes
    ServerAliveInterval 60

Host ada2
    HostName ada-gw2.iiit.ac.in
    User <ada-username>
    IdentityFile ~/.ssh/ada_claude
    IdentitiesOnly yes
    ServerAliveInterval 60
```

Key setup that worked on Windows:

1. `ssh-keygen -t ed25519 -f ~/.ssh/ada_claude` (with a passphrase).
2. Copy the `.pub` with `scp` (password), then append it to `~/.ssh/authorized_keys` on Ada with LF line
   endings: Windows-created `.pub` files have CRLF, and the stray CR made sshd ignore the key until removed.
   `chmod 700 ~/.ssh; chmod 600 ~/.ssh/authorized_keys`.
3. Enable the Windows ssh-agent once (admin PowerShell: `Set-Service ssh-agent -StartupType Automatic;
   Start-Service ssh-agent`), then `ssh-add ~/.ssh/ada_claude` (passphrase once; survives reboots).
4. Use **Windows** `ssh.exe` / `scp.exe` (`C:\Windows\System32\OpenSSH`); Git Bash's own ssh cannot see that agent.
   From Git Bash: `/c/Windows/System32/OpenSSH/ssh.exe ada hostname`.

**Gotchas we hit:**

- Several failed or aborted logins in a row make a gateway refuse connections for ~10–15 min
  (`kex_exchange_identification: Connection closed by remote host`). Don't retry in a loop; try the other gateway.
- Password logins sometimes get `Connection reset` right after the password; just retry once.
- PowerShell 5.1 strips double quotes inside arguments passed to `ssh`, so complex remote one-liners break.
  Put remote logic in a script and send it with `cmd /c "ssh ada bash -s < script.sh"`, or use Git Bash.
- `rsync` is not installed on the Windows laptop: use `git` (the repo is public) for code and `scp` for results.

## Storage plan

| Path | Size | Visible from | Purge | Use for |
|---|---|---|---|---|
| `/home2/$USER` | 50 GB quota (shared with your other projects) | everywhere | never (daily backup) | code (`~/project`), Python env, small results |
| `/share1/$USER` | 100 GB (per guide) | **gateway only** | never | long-term data: TRACE features, our frozen-slide bags |
| `/scratch` | 1.8 TB per node | **that compute node only** | files older than 7 days | per-job working copy, raw slide batches |
| `/ssd_scratch` | 880 GB per node | that compute node only | 7 days | fast I/O for training |

- Compute nodes have internet, so jobs can download **directly** from GDC / Hugging Face to their node's `/scratch`.
- Data on one node's `/scratch` is not visible elsewhere and is purged after 7 days. Anything worth keeping
  goes back to `/share1` (from the gateway) or to home if small.
- Raw frozen slides (430 GB) are downloaded inside the extraction job to `/scratch`, encoded, and deleted.
  Only the resulting bags (~12 GB fp16) are kept.
- TRACE UNI v1 20× features (62 GB): either keep in `/share1` and copy to `/scratch` at job start, or re-download
  from Hugging Face inside the job. Still to test: whether a compute node can `scp` from the gateway.

## Jobs

- Never run heavy work on the gateway.
- Template header (see `jobs/check_env.sbatch`):
  ```
  #SBATCH -A research
  #SBATCH --qos=medium
  #SBATCH -p u22
  #SBATCH --constraint=2080ti
  #SBATCH --gres=gpu:1
  #SBATCH -c 6
  #SBATCH --mem-per-cpu=3G
  #SBATCH --time=04:00:00
  #SBATCH --output=logs/%j.out
  ```
  Default time is 1 h if `--time` is omitted. Max memory 3 GB per CPU.
- **Quota is per user across all projects.** If your other jobs already use most of the 40 CPUs / 4 GPUs,
  new jobs wait with reason `QOSMaxCpuPerUserLimit`. Request fewer CPUs, or a pending job can be shrunk with
  `scontrol update JobId=<id> NumCPUs=<n> MinCPUsNode=<n>`.
- Interactive (debugging): `srun --pty -p u22 -A research --gres=gpu:1 -c 4 --mem-per-cpu=3G --time=02:00:00 bash -l`
- Modules: `module load u22/cuda/12.1` (PyTorch cu121 wheels) and create a venv in home with `u22/python/3.12.4`.
- Monitor: `squeue -u $USER`, `sacct -j <id> --format=JobID,State,Elapsed,MaxRSS,NodeList`, `tail -n 50 logs/<id>.out`.
- Print with `flush=True` so logs update while the job runs.

## Still to check

- [ ] Can a compute node `scp` from the gateway's `/share1`?
- [ ] Is our `/share1/$USER` quota really 100 GB (`du -sh /share1/$USER`, ask hpc.admin if unsure)?
