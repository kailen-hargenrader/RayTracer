@echo off
setlocal
set EXE=build\bin\rt_render.exe
set SCENE=ASCII\Test1.json
set OUT=Output\Test1.ppm

if not exist "%EXE%" (
  echo Executable not found: %EXE%
  echo Run: nmake all
  exit /b 1
)

"%EXE%" "%SCENE%" "%OUT%" --bvh --spp 1 --maxDepth 3 --minThroughput 0.02 --roughSamples 1
if errorlevel 1 exit /b 1

python Code\tools\view_ppm.py "%OUT%"

endlocal


