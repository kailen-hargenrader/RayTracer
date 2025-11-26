@echo off
setlocal
set EXE=build\bin\rt_render.exe
set SCENE=ASCII\Test2.json
set OUT=Output\Test2.ppm

if not exist "%EXE%" (
  echo Executable not found: %EXE%
  echo Run: nmake all
  exit /b 1
)

"%EXE%" "%SCENE%" "%OUT%" --bvh --spp 16 --maxDepth 5 --minThroughput 0.02
if errorlevel 1 exit /b 1

python Code\tools\view_ppm.py "%OUT%"

endlocal