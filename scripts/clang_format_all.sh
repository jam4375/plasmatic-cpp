#!/usr/bin/env bash

git ls-files '*.cpp' | xargs clang-format -style=file -i
git ls-files '*.h' | xargs clang-format -style=file -i