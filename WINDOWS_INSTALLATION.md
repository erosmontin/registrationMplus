# registrationMplus - Windows Build Guide

## Prerequisites

| Component | Version | Download |
|-----------|---------|----------|
| CMake | 4.x+ | https://cmake.org/download/ |
| Visual Studio 2026 Insiders (or 2022) | With C++ workload | https://visualstudio.microsoft.com/downloads/ |
| Boost | 1.87.0+ | https://archives.boost.io/release/1.87.0/source/boost_1_87_0.zip |
| ITK | 4.13 | https://github.com/InsightSoftwareConsortium/ITK/tree/v4.13.0 |
| Git | Latest | https://git-scm.com/download/win |

> **Important:** When installing Visual Studio, make sure to select the **"Desktop development with C++"** workload.

---

## Build Pipeline

### 1. Set Up Directory Structure

```powershell
mkdir C:\Users\%USERNAME%\SRC
mkdir C:\Users\%USERNAME%\BUILD
mkdir C:\Users\%USERNAME%\BUILD\ITK4_13
mkdir C:\Users\%USERNAME%\BUILD\mplus
```

### 2. Download Sources

```powershell
cd C:\Users\%USERNAME%\SRC

# Clone ITK 4.13
git clone -b v4.13.0 https://github.com/InsightSoftwareConsortium/ITK.git

# Download and extract Boost
Invoke-WebRequest -Uri "https://archives.boost.io/release/1.87.0/source/boost_1_87_0.zip" -OutFile boost_1_87_0.zip
Expand-Archive boost_1_87_0.zip -DestinationPath .
```

### 3. Patch ITK 4.13 Source

Open `SRC\ITK\Modules\Core\Common\include\itkImageAlgorithm.h` and find the `StaticCast` struct (around line 190-195).

**Change:**
```cpp
TOutputPixelType operator() (const TInputPixelType i) const
{ return static_cast<TOutputPixelType>(i); }
```

**To:**
```cpp
TOutputPixelType operator() (const TInputPixelType i) const
{ return TOutputPixelType(i); }
```

### 4. Build ITK 4.13

```powershell
cd C:\Users\%USERNAME%\BUILD\ITK4_13
cmake "-G" "Visual Studio 18 2026" "-A" "x64" ..\..\SRC\ITK\ "-DCMAKE_POLICY_VERSION_MINIMUM=3.5" "-DBUILD_TESTING=OFF" "-DBUILD_EXAMPLES=OFF" "-DITKV3_COMPATIBILITY=ON" "-DMODULE_ITKV3Compatibility=ON"
cmake --build . --config Release -- /maxcpucount
```

> If using Visual Studio 2022 instead, replace `"Visual Studio 18 2026"` with `"Visual Studio 17 2022"` in all cmake commands.

#### Verify ITK Build

```powershell
dir C:\Users\%USERNAME%\BUILD\ITK4_13\lib\Release\*.lib | Measure-Object
```

You should see **88+ .lib files**.

### 5. Build Boost (program_options)

Boost must be built from a **Developer Command Prompt** so the compiler is available.

Open **CMD** (not PowerShell):

```cmd
call "C:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Auxiliary\Build\vcvarsall.bat" amd64
cd C:\Users\%USERNAME%\SRC\boost_1_87_0
bootstrap.bat vc143
b2.exe --with-program_options address-model=64 variant=release link=static threading=multi toolset=msvc
```

> For Visual Studio 2022, the vcvarsall.bat path is:
> `"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat"`

### 6. Build registrationMplus

```powershell
cd C:\Users\%USERNAME%\BUILD\mplus
cmake "-G" "Visual Studio 18 2026" "-A" "x64" ..\..\SRC\registrationMplus\src\ "-DCMAKE_POLICY_VERSION_MINIMUM=3.5" "-DITK_DIR=C:\Users\%USERNAME%\BUILD\ITK4_13" "-DBOOST_ROOT=C:\Users\%USERNAME%\SRC\boost_1_87_0"
cmake --build . --config Release -- /maxcpucount
```

---

## Troubleshooting

### CMake: `Compatibility with CMake < 3.5 has been removed`

Add the policy flag to your cmake command:
```
"-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
```

### CMake: `Invalid CMAKE_POLICY_VERSION_MINIMUM value "3"`

PowerShell splits `-D` arguments at the `.` character. **Wrap all `-D` flags in double quotes:**
```powershell
# Wrong:
cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5

# Correct:
cmake "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
```

Alternatively, use **CMD** instead of PowerShell to avoid this issue entirely.

### CMake: `Could not find any instance of Visual Studio`

1. Ensure Visual Studio is installed with the **"Desktop development with C++"** workload
2. Verify installation:
   ```powershell
   & "C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe" -all -property installationPath
   ```
3. If empty, open **Visual Studio Installer** → **Modify** → check **"Desktop development with C++"**

### CMake: `Does not match the generator used previously`

Delete the cache and reconfigure:
```powershell
Remove-Item -Recurse -Force CMakeCache.txt, CMakeFiles
```

### CMake: `Could NOT find Boost`

Ensure `-DBOOST_ROOT` points to your Boost source directory:
```
"-DBOOST_ROOT=C:\Users\%USERNAME%\SRC\boost_1_87_0"
```

If CMake 4.x still can't find it, add these additional flags:
```
"-DBoost_NO_BOOST_CMAKE=ON" "-DCMAKE_POLICY_DEFAULT_CMP0167=OLD"
```

### Boost: `Unknown toolset: vcunk`

Boost doesn't recognize VS 2026. Specify the toolset manually:
```cmd
bootstrap.bat vc143
```

### Boost: `'cl' is not recognized`

You need to run Boost build from a **Developer Command Prompt**, not regular PowerShell. Open CMD and run:
```cmd
call "C:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Auxiliary\Build\vcvarsall.bat" amd64
```

### ITK Build: `error C2143: syntax error in itkImageAlgorithm.h`

The VS 2026 compiler is stricter with `static_cast` on certain types. Apply the patch described in **Step 3** above.

### registrationMplus: `Cannot open include file: 'itkTransformToDeformationFieldSource.h'`

ITK 4.13 moved this header to the V3Compatibility module. Reconfigure and rebuild ITK with:
```
"-DITKV3_COMPATIBILITY=ON" "-DMODULE_ITKV3Compatibility=ON"
```

### General: `std::complex` deprecation warnings (C4996 / STL4037)

These are non-fatal warnings. To silence them, add to your cmake configure:
```
"-DCMAKE_CXX_FLAGS=/D_SILENCE_NONFLOATING_COMPLEX_DEPRECATION_WARNING"
```

---

## Visual Studio Version Reference

| Visual Studio Version | Generator Name | vcvarsall.bat Path |
|----------------------|----------------|--------------------|
| VS 2026 Insiders | `Visual Studio 18 2026` | `C:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Auxiliary\Build\vcvarsall.bat` |
| VS 2022 | `Visual Studio 17 2022` | `C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat` |
| VS 2019 | `Visual Studio 16 2019` | `C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat` |

---

## Notes

- All deprecation warnings from CMake (`Compatibility with CMake < 3.10 will be removed...`) are **safe to ignore**
- Warning flags showing `Failed` during ITK configuration (e.g., `-Wno-uninitialized - Failed`) are normal — those are GCC/Clang flags that don't apply to MSVC
- Build times: ITK ~20-60 min, Boost ~2-5 min, registrationMplus ~5-10 min