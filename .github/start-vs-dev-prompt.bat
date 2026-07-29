@echo off

set "ARCH=%~1"

if "%ARCH%"=="" (
  set "ARCH=x64"
)

rem https://github.com/microsoft/vswhere/wiki/Start-Developer-Command-Prompt
for /f "usebackq delims=" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -prerelease -latest -property installationPath`) do (
  if exist "%%i\Common7\Tools\vsdevcmd.bat" (
    call  "%%i\Common7\Tools\vsdevcmd.bat" -arch=%ARCH%
    exit /b
  )
)

echo Instance or command prompt not found
exit /b 2
