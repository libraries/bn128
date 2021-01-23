BIN_PREFIX=/root/app/riscv/bin/riscv64-unknown-elf
GPP=$BIN_PREFIX-g++

set -ex

mkdir -p build

clang-format -i -style="{ColumnLimit: 120}" include/*
clang-format -i -style="{ColumnLimit: 120}" test/*

g++ -fno-exceptions -Os -Iinclude -I/src/intx/include -o build/test test/test_bn128_mont.cpp
./build/test
echo "ok"

$GPP -fno-exceptions -Os -march=rv64gc -Iinclude -I/src/intx/include -o build/test test/test_bn128_mont.cpp
pm c /src/ckb_vm_run/target/release/int64 build/test
pm l int64
