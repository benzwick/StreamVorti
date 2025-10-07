# Installing StreamVorti on Setonix with Spack

Building StreamVorti on [Setonix](https://pawsey.org.au/systems/setonix/) using your own Spack (newer than Pawsey's 0.21.0).

## Key Setonix Requirements

- **Project**: `pawsey1243`
- **Build location**: `/software/projects/pawsey1243/shared/` (shared by all project members) - [docs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925878/How+to+Manually+Build+Software)
- **Use compute nodes**: NOT login nodes - [docs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925954/Compiling)
- **Use `sg`**: All `/software` operations need `sg pawsey1243 -c 'command'` - [docs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925886/Spack)
- **Partition**: Use `work` (NOT `debug`) for builds - [docs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix)

## Shared Directory Permissions

The `shared` directory setup ensures all project members can collaborate:

1. Create directory with `sg pawsey1243 -c 'mkdir -p shared'`
2. Set permissions with `chmod 2775` which gives:
   - `2` = setgid bit (inherits group)
   - `7` = owner rwx
   - `7` = group rwx (read+write+execute)
   - `5` = others rx (read+execute)
3. This results in `drwxrwsr-x` where group members can read and write

All new files/folders created inside will automatically get the `pawsey1243` group and be accessible to all project members.

## Installation Steps

```bash
# 1. Create shared directory with correct permissions
# See: https://pawsey.atlassian.net/wiki/spaces/US/pages/51925886/Spack (sg usage)
# See: https://pawsey.atlassian.net/wiki/spaces/US/pages/51925878/How+to+Manually+Build+Software (build location)
cd /software/projects/pawsey1243/
sg pawsey1243 -c 'mkdir -p shared'
sg pawsey1243 -c 'chmod 2775 shared'
# Verify: ls -la should show drwxrwsr-x with pawsey1243 group

# 2. Request compute node (see: https://pawsey.atlassian.net/wiki/spaces/US/pages/51925964/Job+Scheduling)
# Example from Pawsey docs for shared node access:
salloc -p work -N 1 -n 1 -c 64 --mem=115G --time=4:00:00 -A pawsey1243
# Or simpler for full node:
# salloc --nodes=1 --partition=work --time=4:00:00 --account=pawsey1243
#
# When allocation succeeds, you'll see:
#   salloc: Granted job allocation XXXXXXX
#   salloc: Waiting for resource configuration
#   salloc: Nodes nidXXXXXX are ready for job
# Your shell prompt will change from setonix-XX to nidXXXXXX showing you're on the compute node.
# All subsequent commands now execute on the compute node (not login node).
# Use 'exit' to release the node when done.

# 3. Clone StreamVorti to shared directory
# See: https://pawsey.atlassian.net/wiki/spaces/US/pages/51925886/Spack (sg for /software operations)
cd /software/projects/pawsey1243/shared/
sg pawsey1243 -c 'git clone https://github.com/benzwick/StreamVorti.git streamvorti'
cd streamvorti

# 4. Install Spack with Pawsey's Setonix config
# See: https://github.com/PawseySC/pawsey-spack-config (Pawsey's Spack configs)
# See: https://spack.readthedocs.io/en/latest/getting_started.html (Spack installation)
sg pawsey1243 -c './scripts/build-with-spack/02-install-spack.sh ../spack'
sg pawsey1243 -c 'git clone https://github.com/PawseySC/pawsey-spack-config.git ../pawsey-spack-config'
# Copy compilers.yaml (optimized flags and Cray compiler paths)
cp ../pawsey-spack-config/systems/setonix/configs/site/compilers.yaml ../spack/etc/spack/
# Copy minimal packages.yaml (compiler preferences and providers)
# Note: This is a simplified version of Pawsey's packages.yaml without version constraints
# See scripts/build-with-spack/setonix-packages.yaml for details
cp scripts/build-with-spack/setonix-packages.yaml ../spack/etc/spack/packages.yaml
# Configure Spack to use /software for caching (not home directory which has 1GB limit)
export SPACK_USER_CACHE_PATH=/software/projects/pawsey1243/shared/spack-cache
mkdir -p $SPACK_USER_CACHE_PATH

# 5. Configure compilers
# See: https://spack.readthedocs.io/en/latest/getting_started.html#compiler-configuration
# See: https://pawsey.atlassian.net/wiki/spaces/US/pages/51929054/Setonix+Software+Environment (Cray compilers)
../spack/bin/spack compiler find

# 6. Build StreamVorti
# See: https://spack.readthedocs.io/en/latest/environments.html (Spack environments)
# See: https://pawsey.atlassian.net/wiki/spaces/US/pages/51925954/Compiling (compiling on compute nodes)
sg pawsey1243 -c './scripts/build-with-spack/04-create-spack-environment.sh ../spack streamvorti spack.yaml'
sg pawsey1243 -c './scripts/build-with-spack/07-build-streamvorti.sh ../spack streamvorti'

# 7. Test
./scripts/build-with-spack/08-run-ctest.sh ../spack streamvorti
./scripts/build-with-spack/09-check-built-artifacts.sh
./scripts/build-with-spack/10-run-example-tests.sh ../spack streamvorti
```

## Spack Configuration

Use [Pawsey's Spack config for Setonix](https://github.com/PawseySC/pawsey-spack-config/tree/master/systems/setonix/configs) which includes:
- Cray/GNU/AMD compiler setup
- System MPI and performance libraries
- Optimized directory structure
- Module generation settings

See [Spack configuration docs](https://spack.readthedocs.io/en/latest/configuration.html) for customization.

## Troubleshooting

### Disk quota exceeded errors

If you get "Disk quota exceeded" errors, files may be created with wrong group ownership. Check quota:

```bash
lfs quota -h /software/projects/pawsey1243
```

If your user quota is high (>1GB) instead of project quota, files have wrong group. Fix by:

```bash
# Delete everything and start over
cd /software/projects/pawsey1243/
sg pawsey1243 -c 'rm -rf shared'

# Recreate with correct permissions from step 1
```

### Remove Spack environment and start over

```bash
sg pawsey1243 -c '../spack/bin/spack env rm -y streamvorti'
```

Then re-run step 6 from Installation Steps.

### Clean up failed package repo clones

```bash
rm -rf ~/.spack/package_repos
```

Ensure `SPACK_USER_CACHE_PATH` is set (see step 4).

### Check available package versions

```bash
../spack/bin/spack versions <package-name>
```

## Documentation Links

- [Setonix User Guide](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925434/Setonix+User+Guide)
- [Spack at Pawsey](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925886/Spack)
- [How to Manually Build Software](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925878/How+to+Manually+Build+Software)
- [Compiling on Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925954/Compiling)
- [Filesystem Policies](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925880/Filesystem+Policies)
- [Running Jobs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix)
