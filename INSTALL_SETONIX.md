# Installing StreamVorti on Setonix with Spack

Building StreamVorti on [Setonix](https://pawsey.org.au/systems/setonix/) using your own Spack (newer than Pawsey's 0.21.0).

## Key Setonix Requirements

- **Build location**: `/software/projects/<project>/<user>/` ([docs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925878/How+to+Manually+Build+Software))
- **Use compute nodes**: NOT login nodes ([docs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925954/Compiling))
- **Use `sg`**: All `/software` operations need `sg <project> -c 'command'` ([docs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925886/Spack))
- **Partition**: Use `work` (NOT `debug`) for builds ([docs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix))

## Installation Steps

```bash
# 1. Request compute node
salloc --nodes=1 --partition=work --time=4:00:00 --account=<project>

# 2. Clone StreamVorti to software directory
cd /software/projects/<project>/<user>/
sg <project> -c 'git clone https://github.com/benzwick/StreamVorti.git streamvorti'
cd streamvorti

# 3. Install Spack with Pawsey's Setonix config
sg <project> -c './scripts/build-with-spack/02-install-spack.sh ../spack'
sg <project> -c 'git clone https://github.com/PawseySC/pawsey-spack-config.git ../pawsey-spack-config'
cp -r ../pawsey-spack-config/systems/setonix/configs/site/* ../spack/etc/spack/

# 4. Configure compilers
../spack/bin/spack compiler find

# 5. Build StreamVorti
sg <project> -c './scripts/build-with-spack/04-create-spack-environment.sh ../spack streamvorti spack.yaml'
sg <project> -c './scripts/build-with-spack/07-build-streamvorti.sh ../spack streamvorti'

# 6. Test
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

## Documentation Links

- [Setonix User Guide](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925434/Setonix+User+Guide)
- [Spack at Pawsey](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925886/Spack)
- [How to Manually Build Software](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925878/How+to+Manually+Build+Software)
- [Compiling on Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925954/Compiling)
- [Filesystem Policies](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925880/Filesystem+Policies)
- [Running Jobs](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix)
