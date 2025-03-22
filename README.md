# ODE

## Build

```bash
cmake -B <tmp-build-directory> -G "Ninja" -D CMAKE_BUILD_TYPE=Release -S <project-directory>
cmake --build <tmp-build-directory>
```

Optionnal:

- `-D CMAKE_C_COMPILER=gcc`
- `-D CMAKE_CXX_COMPILER=g++`
