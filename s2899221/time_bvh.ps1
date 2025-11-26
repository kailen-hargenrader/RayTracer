param()

$ErrorActionPreference = 'Stop'

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root

$exe = Join-Path $root 'build\bin\rt_render.exe'
$scene = Join-Path $root 'ASCII\Test2.json'
$outNo = Join-Path $root 'Output\Test2_nobvh.ppm'
$outBvh = Join-Path $root 'Output\Test2_bvh.ppm'

$spp = 4
$depth = 3
$rough = 1

if (-not (Test-Path $exe)) {
  Write-Host "Executable not found: $exe"
  Write-Host "Build the project first (e.g., nmake all)"
  exit 1
}

if (-not (Test-Path (Join-Path $root 'Output'))) {
  New-Item -ItemType Directory -Path (Join-Path $root 'Output') | Out-Null
}

Write-Host "Benchmarking `"$scene`" with SPP=$spp maxDepth=$depth roughSamples=$rough"

$msNo = (Measure-Command { & $exe $scene $outNo --spp $spp --maxDepth $depth --roughSamples $rough | Out-Null }).TotalMilliseconds
$msBvh = (Measure-Command { & $exe $scene $outBvh --bvh --spp $spp --maxDepth $depth --roughSamples $rough | Out-Null }).TotalMilliseconds

$sNo = [math]::Round($msNo/1000.0, 2)
$sBvh = [math]::Round($msBvh/1000.0, 2)
$speedup = "{0:N2}" -f ($msNo / $msBvh)

Write-Host ""
Write-Host "Results:"
Write-Host ("  No BVH : {0} ms  (~{1}s)" -f $msNo, $sNo)
Write-Host ("  With BVH: {0} ms  (~{1}s)" -f $msBvh, $sBvh)
Write-Host ("  Speedup : {0}x" -f $speedup)


