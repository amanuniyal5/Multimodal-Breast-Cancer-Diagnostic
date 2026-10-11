# Working on Ada (IIIT HPC)

Summary of the official [Ada User Guide](https://hpc.iiit.ac.in/wiki/index.php?title=Ada_User_Guide) plus what we
learned connecting to it. Needs the IIIT VPN off campus.

## Connecting (changed in 2026)

`ada.iiit.ac.in` no longer accepts SSH. Use the gateways **`ada-gw1.iiit.ac.in`** (10.4.24.30) or
**`ada-gw2.iiit.ac.in`** (10.4.24.29). Password login works once; after that use a key.

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

Install the key once (Windows has no `ssh-copy-id`): copy the `.pub` file with `scp`, then log in with your
password and append it to `~/.ssh/authorized_keys` (`chmod 700 ~/.ssh; chmod 600 ~/.ssh/authorized_keys`).

**Gotchas we hit:**

- Several failed or aborted logins in a row make a gateway refuse connections for a while
  (`kex_exchange_identification: Connection closed by remote host`). Wait ~10–15 min; use the other gateway.
- PowerShell 5.1 strips double quotes inside arguments passed to `ssh`, so complex remote one-liners break.
  Put remote logic in a script file on Ada, or run `ssh` from Git Bash.
- `rsync` is not installed on the Windows laptop: use `git` (the repo is public) for code and `scp` for results.

## Storage (important for our data plan)

| Path | Quota | Visible from | Backup | Purge | Use for |
|---|---|---|---|---|---|
| `/home/$USER` | 25 GB | master + compute | daily | never | code, conda env, small results |
| `/share1/$USER` | 100 GB | **master only** | no | never | long-term data: TRACE features, our frozen-slide bags |
| `/scratch` | none (1.8 TB per node) | **only that compute node** | no | **files older than 7 days** | per-job working copy of data, raw slide batches |
| `/ssd_scratch` | none (960 GB per node) | only that compute node | no | 7 days | fast I/O for training |

Consequences:

- Data stored on one node's `/scratch` is **not** visible to other nodes and disappears after 7 days.
- Long-term data lives in `/share1`. A job copies what it needs to the node at start:
  `scp ada:/share1/$USER/<file> /scratch/$USER/` (the guide's pattern), and copies outputs back at the end.
- Budget for `/share1` (100 GB): TRACE UNI v1 20× features 62 GB (or 10× 16 GB) + frozen-slide bags ~12 GB
  (fp16) + CONCH 8 GB. Raw frozen slides (430 GB) are downloaded inside a job to `/scratch`, encoded,
  then deleted; they never go to `/share1`.

## Partitions, account, QoS

| Partition | Nodes | GPUs | Max time |
|---|---|---|---|
| `short` | gnode01–02 | 4 per node | 6 h (debugging) |
| `long` | gnode03–92 (01–40: GTX 1080 Ti, 43–92: RTX 2080 Ti) | 4 per node | unlimited |

- Account `research` (everyone), default QoS `medium`: up to 40 CPUs, **4 GPUs per user**, 4 running jobs,
  8 submitted, 4-day wall time. QoS `low`: 1 GPU, 1 job.
- Account `sub`, QoS `sub`: idle nodes only, 6 h, low priority, up to 4 GPUs.
- Keep about 1 GPU : 10 CPUs. Default memory 1 GB per CPU, maximum 3 GB per CPU.
- Check your associations: `sacctmgr show assoc user=$USER format=Account,QOS,DefaultQOS`.

## Modules

`module avail`, then for example `module load u18/cuda/11.6 u18/cudnn/8.4.0-cuda-11.6` and
`module load u18/python/3.10.2`. Custom Python environments must live in your home directory.

## Jobs

- Never run heavy work on the gateway/login node.
- Interactive (max 6 h): `srun --pty -p long -A research --gres=gpu:1 -c 10 --mem-per-cpu=2G bash -l`
- Batch: `sbatch jobs/<script>.sbatch` (templates in `jobs/`). Print with `flush=True` so logs update.
- Monitor: `squeue -u $USER`, `sacct -j <id> --format=JobID,State,Elapsed,MaxRSS,NodeList`,
  `tail -n 50 logs/<id>.out`.
- `CUDA_VISIBLE_DEVICES` holds the GPUs assigned to the job.

## Open questions (fill in once checked)

- [ ] Home path on the gateways (`/home/$USER` or `/home2/$USER`) and whether both gateways share it.
- [ ] Do compute nodes have internet access (needed to download from GDC / Hugging Face inside jobs)?
- [ ] How to request RTX 2080 Ti nodes specifically (`sinfo -o "%N %f %G"`: node features / GRES types).
