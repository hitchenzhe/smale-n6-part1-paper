param([int]$MemoryMB=1024)
$ErrorActionPreference='Stop'
$reviewDir=$PSScriptRoot
$verificationDir=Split-Path -Parent $reviewDir
$pythonExe=Join-Path $verificationDir '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $pythonExe)) { throw 'Run setup_python.ps1 next to this wrapper first.' }
$jobs=@(
 (Join-Path $verificationDir 'verify_n6_structure.py'),
 'generate_no8_revised_order_certificate.py',
 'generate_shared_energy_certificate.py',
 'verify_no20_energy.py',
 'generate_no50_resonance_certificate.py',
 'generate_no79_bridge_certificate.py',
 'verify_no79_shape_energy.py',
 'generate_no79_sturm_certificate.py',
 'verify_no79_hurwitz.py',
 'generate_no52_certificate.py',
 'verify_no52_valuation_partition.py',
 'verify_revised_core_claims.py',
 'verify_no66_function_field_cover.py',
 'generate_no88_collinear_certificate.py',
 'generate_no88_fiber_certificate.py',
 'generate_no88_vertical_certificate.py',
 'generate_no88_smooth_certificate.py',
 'generate_no88_zero_star_certificate.py',
 'verify_no88_root_norm_identity.py',
 'verify_no88_CM_edge_quadratic.py',
 'generate_crossed_torque_certificate.py',
 'generate_no93_interval_certificate.py',
 'generate_no98_interval_certificate.py',
 'generate_no63_certificate.py',
 'derive_no94_finite_fiber.py',
 'generate_no94_structure_certificate.py',
 'generate_no94_equal_spoke_certificate.py',
 'generate_no94_projective_spoke_cover.py',
 'generate_no94_local_certificates.py',
 'verify_no94_distinct_spoke_fiber.py',
 'verify_no94_h3_fiber.py'
)
foreach($job in $jobs) {
 if (-not [IO.Path]::IsPathRooted($job)) { $job=Join-Path $reviewDir $job }
 Write-Host ('Replaying ' + [IO.Path]::GetFileName($job))
 & $pythonExe (Join-Path $reviewDir 'run_certificate_job.py') $job --memory-mb $MemoryMB --timeout-s 900
 if($LASTEXITCODE -ne 0) { throw ('Replay failed: ' + $job) }
}
& $pythonExe (Join-Path $reviewDir 'summarize_appendix_replay.py')
if($LASTEXITCODE -ne 0) { throw 'Appendix replay summary contains failed or missing jobs.' }
Write-Host 'Current implemented checks passed. Consult appendix_replay_manifest.json for partial and missing proof obligations.'
