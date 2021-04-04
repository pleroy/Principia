# Instructions for building Principia on Windows

Before starting, make sure the following are installed on your machine:
* Visual Studio 2019 version 16.8.1 (note that 16.9.0 does not appear to work) with C# and C++ support;
* Windows SDK 10.0.18362;
* .NET Framework 4.7.2 SDK (and Targeting Pack).

The solution contains a C# project named `coverage_analyser` which requires
the Enterprise edition of Visual Studio to build.  It's only a development
tool for our own use, so if you use a different edition, just do 
`Project > Unload Project` on that project.  You will still be able to build 
and test the mod.

These instructions use the git Powershell provided by [GitHub for Windows](https://windows.github.com/).
We assume a working installation of Kerbal Space Program version 1.8.1 is found in `<KSP directory>`.

The repository is found at https://github.com/mockingbirdnest/Principia.git.
Pick a directory `<root>` in which you will install Principia and its
dependencies.
This directory should not contain any of the following subfolders:
- `Principia`;
- `KSP Assemblies`;
- `Google`.

This project depends upon:
- the KSP assembly `Assembly-CSharp.dll`, found in `<KSP directory>\KSP_x64_Data\Managed`;
- the Unity assemblies `UnityEngine.CoreModule.dll`, `UnityEngine.dll`, `UnityEngine.ImageConversionModule.dll`, `UnityEngine.IMGUIMode.dll`, `UnityEngine.InputLegacyModule.dll`, `UnityEngine.PhysicsModule.dll`, `UnityEngine.TextRenderingModule.dll` and `UnityEngine.UI.dll`, found in
  `<KSP directory>\KSP_x64_Data\Managed`;
- our [fork](https://github.com/mockingbirdnest/glog) of the Google glog
  library;
- our [fork](https://github.com/mockingbirdnest/googletest) of the Google googletest
  and googlemock libraries;
- our [fork](https://github.com/mockingbirdnest/protobuf) of the Google
  protobuf library;
- our [fork](https://github.com/mockingbirdnest/benchmark) of the Google
  benchmark library;
- our [fork](https://github.com/mockingbirdnest/gipfeli) of the Google gipfeli library;
- our [fork](https://github.com/mockingbirdnest/abseil-cpp) of the Google Abseil C++ library;
- our [fork](https://github.com/mockingbirdnest/zfp) of the LLNL zfp library;
- parts of the Chromium codebase (for stack tracing support in glog on Windows),
  *modified according to the instructions below*.

## Installation steps

### Dowloading Principia

In `<root>`, run `git clone https://github.com/mockingbirdnest/Principia.git`.

### KSP and Unity assemblies

In order to build for KSP 1.8.1, copy the corresponding KSP 1.8.1 assemblies to `<root>\KSP Assemblies\1.8.1`

### Downloading the dependencies

In `<root>`, run the following commands.
```powershell
mkdir "Google"
mkdir "Third Party"
push-location -path "Google"
git clone "https://chromium.googlesource.com/chromium/src.git" chromium -n --depth 1 -b "40.0.2193.1"
$GitPromptSettings.RepositoriesInWhichToDisableFileStatus += join-path  (gi -path .).FullName chromium
push-location -path "chromium"
git config core.sparsecheckout true
copy "..\..\Principia\documentation\setup files\chromium_sparse_checkout.txt" ".git/info/sparse-checkout"
git checkout
copy "..\..\Principia\documentation\setup files\chromium.patch"
git am "chromium.patch"
rm "chromium.patch"
pop-location
pop-location
```
### Building

In `<root>`, run the following command:
```powershell
.\Principia\rebuild_all_solutions.ps1
```

# Instructions for building Principia on Linux and macOS

*Note that the released binaries for Linux and macOS are built using [Azure pipelines](https://dev.azure.com/mockingbirdnest/Principia/_build).  The instructions below are best effort.*

## Build Prerequisites
Before starting, make sure the following are installed on your machine:
### Linux
* Plugin build prerequisites: `build-essential`, `clang`, `libc++-dev`, `libc++abi-dev`, `subversion`, and `git`;
* Adapter build prequisites: `msbuild`
* Runtime dependencies: `libc++1`.
### macOS
* You will need `clang` version 8 and the Xcode command-line tools.
Clang 8 is included as part of Xcode versions 11.0 through 11.3.1.
If you can't install Xcode 11 (possibly because you are using a macOS version prior to 10.14.4),
you can install clang manually using [MacPorts](https://www.macports.org) or [Homebrew](https://brew.sh).
You will still need to install some version of Xcode for its command-line tools.
* To build the adapter, you will need `msbuild`. It is included as part of [Mono](https://www.mono-project.com/download/stable/), and can alternatively be installed directly using a package manager.
You might be able to build it using [Visual Studio for Mac](https://visualstudio.microsoft.com/vs/mac/),
but this is untested.

## Installing the dependencies

In `<root>/Principia`, run the following command:
```bash
./install_deps.sh
```
This will install and compile all the third-party components that Principia uses.  Don't proceed with the next step unless this step has completed without errors.

### KSP and Unity assemblies
In order to build the adapter, you will need references to the KSP and Unity assemblies.
Principia expects to find these in a directory adjacent to the principia directory.
For example, in order to build for KSP 1.10.1, the assembly directory should be linked to `../KSP Assemblies/1.10.1` (relative to the Principia directory).

A shell script is provided that will do this automatically on macOS:
```
./add_ksp_assemblies_macos.sh
```

## Building Principia

In `<root>/Principia`, run the following command:
```bash
make
```
See the Makefile for more options.
If some of the unit tests fail, you may or may not be able to run the resulting version of Principia.
