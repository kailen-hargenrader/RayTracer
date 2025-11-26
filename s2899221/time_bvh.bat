@echo off
setlocal

rem Ensure we run relative paths from this script's directory and delegate to PowerShell script for reliability
pushd "%~dp0"

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0time_bvh.ps1"
set EXITCODE=%ERRORLEVEL%

popd
endlocal
exit /b %EXITCODE%


