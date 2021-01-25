set -ex

GPP=g++
GPP_RISCV=$RISCV/bin/riscv64-unknown-elf-g++

mkdir -p build

$GPP -fno-exceptions -Os -Iinclude -Iintx/include -o build/test test/test_bn128.cpp
./build/test
echo "ok"

$GPP_RISCV -fno-exceptions -Os -march=rv64gc -Iinclude -Iintx/include -o build/test test/test_bn128.cpp
./ckb-vm-run/target/release/asm build/test
