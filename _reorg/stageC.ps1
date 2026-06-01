$ErrorActionPreference = 'Stop'
function Md5($p){ (Get-FileHash -Algorithm MD5 -LiteralPath $p).Hash }

$libHashes = Get-Content '_reorg\lib_hashes.json' -Raw | ConvertFrom-Json
$log = New-Object System.Collections.ArrayList
@(Get-Content '_reorg\reorg_actions.json' -Raw | ConvertFrom-Json) | ForEach-Object { [void]$log.Add($_) }

function CopyFromLib($files,$dst,$group){
  if(-not (Test-Path -LiteralPath $dst)){ New-Item -ItemType Directory -Force -Path $dst | Out-Null }
  foreach($f in $files){
    $s = Join-Path 'lib' $f
    $d = Join-Path $dst $f
    if(Test-Path -LiteralPath $d){ throw "UNEXPECTED existing dest: $d" }
    $b = Md5 $s
    Copy-Item -LiteralPath $s -Destination $d
    $a = Md5 $d
    if($a -ne $b -or $a -ne $libHashes.$f){ throw "copy-from-lib mismatch: $f -> $dst" }
    [void]$log.Add([ordered]@{ group=$group; action='copy-from-lib'; source=$s; dest=$d; before=$b; after=$a })
  }
  Write-Host ("LIB->{0}  [{1} files]" -f $dst,$files.Count)
}

function CopyAsIs($files,$srcDir,$dst,$group){
  if(-not (Test-Path -LiteralPath $dst)){ New-Item -ItemType Directory -Force -Path $dst | Out-Null }
  foreach($f in $files){
    $s = Join-Path $srcDir $f
    $d = Join-Path $dst $f
    if(Test-Path -LiteralPath $d){ throw "UNEXPECTED existing dest: $d" }
    $b = Md5 $s
    Copy-Item -LiteralPath $s -Destination $d
    $a = Md5 $d
    if($a -ne $b){ throw "copy-as-is mismatch: $f" }
    [void]$log.Add([ordered]@{ group=$group; action='copy-as-is'; source=$s; dest=$d; before=$b; after=$a })
  }
  Write-Host ("ASIS->{0}  [{1} files]" -f $dst,$files.Count)
}

# 04_interaction/pipeline : canonical interaction stack (10 from lib)
$pipe = 'ls_basis.R','spatial_utils.R','marginal_utils.R','gibbs_stage_c_full.R','fit_spatial_reml.R','canonical_schema.R','ls_interaction.R','ls_interaction_core.cpp','gibbs_interaction.R','fit_ours_interaction_wrapper.R'
CopyFromLib $pipe '04_interaction\pipeline' '04_interaction/pipeline/'

# 04_interaction/scenarioA_compare : 6 shared from lib + 11 unique copied as-is
$saShared = 'canonical_schema.R','fit_spatial_reml.R','gibbs_stage_c_full.R','ls_basis.R','marginal_utils.R','spatial_utils.R'
CopyFromLib $saShared '04_interaction\scenarioA_compare' '04_interaction/scenarioA_compare/'
$saSrc = 'archive\interaction_old\comparison_stage1\scenario_A'
$saUnique = 'compute_recovery_metrics.R','fit_inla_wrapper.R','fit_mgcv_wrapper.R','fit_ours_wrapper.R','run_stage1_scenarioA_multiseed.R','scenarioA_settings.md','simulate_scenario_A.R','test_fit_inla.R','test_fit_mgcv.R','test_fit_ours.R','test_simulate_scenario_A.R'
CopyAsIs $saUnique $saSrc '04_interaction\scenarioA_compare' '04_interaction/scenarioA_compare/'

# Inject missing lib working copies into live experiment folders
CopyFromLib @('ls_basis.R','spatial_utils.R') '02_reml' '02_reml/'
CopyFromLib @('ls_basis.R','spatial_utils.R','fit_spatial_reml.R','marginal_utils.R','baby_bayes.R','gibbs_bayes.R') '03_bayes_main\n400' '03_bayes_main/n400/'
CopyFromLib @('ls_basis.R','spatial_utils.R','fit_spatial_reml.R','marginal_utils.R','baby_bayes.R','gibbs_bayes.R','gibbs_stage_c_full.R') '03_bayes_main\n1000' '03_bayes_main/n1000/'
CopyFromLib @('ls_basis.R','spatial_utils.R','fit_spatial_reml.R','marginal_utils.R','baby_bayes.R','gibbs_bayes.R','gibbs_stage_c_full.R') '05_munich_rent\hpc\codes' '05_munich_rent/hpc/codes/'

($log | ConvertTo-Json -Depth 5) | Out-File -Encoding utf8 '_reorg\reorg_actions.json'
"stageC done. total log rows: $($log.Count)"
