# Pathtracer

## Linux/macOS Instructions

### Installation

1. Navigate to the [Embree download page](https://www.embree.org/downloads.html) and follow the instructions under **Linux RPMs** or **macOS PKG Installer** to install Embree
    - As part of this process, you will need to install Intel's Threading Building Blocks (TBB) library. If you do this manually, make sure to set the correct environment variables (`LD_LIBRARY_PATH` for Linux and `DYLD_LIBRARY_PATH` for macOS). You can do this by `source`ing the `tbbvars.[c]sh` shell script in your shell startup script (eg. `.bashrc`). See https://software.intel.com/en-us/node/505529.
2. Navigate to the repository: `cd <repository path>`
3. Create a build directory for CMake: `mkdir build && cd build`
4. Initialize CMake: `cmake ..`

### Usage

1. Build/rebuild the project with `make -j`
2. Render a scene with `./pathtracer <scene file>`
3. The final image will be saved to your current directory

## Windows Instructions

### Installation

1. Navigate to the [Embree download page](https://www.embree.org/downloads.html) and follow the instructions under **Windows MSI Installer** to install Embree
2. Open the Environment Variables window (*Control Panel -> System and Security -> System -> Advanced system settings -> Environment Variables...*). Select the *Path* variable and press *Edit...*. Create a new entry with the path `C:\Program Files\Intel\Embree3 x64\bin`.
3. Checkout the repository from GitHub and initialize CMake.
    - **Visual Studio 2019**: Select *Clone or check out code*, then enter the repository URL.
    - **Visual Studio 2017**: Select *File -> Open -> Open from Source Control*. In the Team Explorer panel, select *Clone*, enter the repository URL, and click *Clone*.
4. Under *Select Startup Item...*, select *pathtracer_win.exe*
5. (Optional) Change *x64-Debug* to *x64-Release* for best performance

### Usage

1. Run *pathtracer_win.exe*
2. Select a scene file when prompted to start rendering
3. The final image will be saved to the same directory as the scene file
    - Visual Studio may not display the image file in the Solution Explorer. If this is the case, simply right-click the directory and select *Open Folder in File Explorer* to view the directory.
