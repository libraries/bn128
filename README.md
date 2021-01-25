# ALT BN128

This repo implements the following EIPs on CKB-VM:

- https://github.com/ethereum/EIPs/blob/master/EIPS/eip-196.md
- https://github.com/ethereum/EIPs/blob/master/EIPS/eip-197.md

The way to use it is dead simple. You can compile riscv or native results, it depends on you.

**Build for CKB/RISC-V**

```sh
# Build test
$ riscv64-unknown-elf-g++ -fno-exceptions -Os -march=rv64gc -Iinclude -Iintx/include -o build/test test/test_bn128.cpp

# Run test on CKB-VM
$ cd ckb-vm-run
$ cargo build --release
$ cd ..
$ ./ckb-vm-run/target/release/asm build/test
```

**Build for native**

```sh
# Build test
$ g++ -fno-exceptions -Os -Iinclude -Iintx/include -o build/test test/test_bn128.cpp

# Run test on native
$ ./build/test
# If the return code is 0, it means that the test has all passed.
$ echo $?
```

**Or use the script directly**

```sh
$ export RISCV=/path/to/riscvdir

$ ./script/test.sh
$ ./script/test_benchmark.sh
```

Read the test code to see how to call `add()`, `mul()` and `pairing()`.
