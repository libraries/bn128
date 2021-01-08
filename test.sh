BIN_PREFIX=/root/app/riscv/bin/riscv64-unknown-elf
GPP=$BIN_PREFIX-g++

set -ex

mkdir -p build

g++ -fno-exceptions -Os -Iinclude -I/src/intx/include -o build/test test/test_bn128.cpp
./build/test
echo "ok"

$GPP -g -fno-exceptions -Os -march=rv64gc -Iinclude -I/src/intx/include -o build/test test/test_bn128.cpp
pm c /src/ckb_vm_run/target/debug/int64 build/test
