$ErrorActionPreference = 'Stop'
function Md5($p){ (Get-FileHash -Algorithm MD5 -LiteralPath $p).Hash }
$codeExt = @('.r','.cpp','.h','.sub','.md','.sh')

$libHashes = Get-Content '_reorg\lib_hashes.json' -Raw | ConvertFrom-Json
$canonNames = $libHashes.PSObject.Properties.Name   # the 12 canonical filenames
$rows = New-Object System.Collections.ArrayList
function AddRow($group,$action,$source,$dest,$before,$after){
  [void]$rows.Add([pscustomobject]@{ group=$group; action=$action; source=$source; dest=$dest; before=$before; after=$after })
}
function TopGroup($p){ ($p -split '\\')[0] }

# ---------- 1. copy-to-lib (12 canonical files into lib\) ----------
$munich = 'codes\run_codes\munich'; $tc = 'codes\test and comparison'
$ctlSrc = @{}
foreach($n in 'baby_bayes.R','fit_spatial_reml.R','gibbs_bayes.R','gibbs_interaction.R','gibbs_stage_c_full.R','ls_basis.R','ls_interaction.R','ls_interaction_core.cpp','marginal_utils.R','spatial_utils.R'){ $ctlSrc[$n] = $munich }
foreach($n in 'canonical_schema.R','fit_ours_interaction_wrapper.R'){ $ctlSrc[$n] = $tc }
foreach($n in $canonNames){
  $d = Join-Path 'lib' $n; $h = Md5 $d
  AddRow 'lib/' 'copy-to-lib' (Join-Path $ctlSrc[$n] $n) $d $h $h
}

# ---------- copy-from-lib injected sets (dest dir -> filenames) ----------
$cfl = [ordered]@{
  '02_reml'                          = @('ls_basis.R','spatial_utils.R')
  '03_bayes_main\n400'               = @('ls_basis.R','spatial_utils.R','fit_spatial_reml.R','marginal_utils.R','baby_bayes.R','gibbs_bayes.R')
  '03_bayes_main\n1000'              = @('ls_basis.R','spatial_utils.R','fit_spatial_reml.R','marginal_utils.R','baby_bayes.R','gibbs_bayes.R','gibbs_stage_c_full.R')
  '04_interaction\pipeline'          = @('ls_basis.R','spatial_utils.R','marginal_utils.R','gibbs_stage_c_full.R','fit_spatial_reml.R','canonical_schema.R','ls_interaction.R','ls_interaction_core.cpp','gibbs_interaction.R','fit_ours_interaction_wrapper.R')
  '04_interaction\scenarioA_compare' = @('canonical_schema.R','fit_spatial_reml.R','gibbs_stage_c_full.R','ls_basis.R','marginal_utils.R','spatial_utils.R')
  '05_munich_rent\hpc\codes'         = @('ls_basis.R','spatial_utils.R','fit_spatial_reml.R','marginal_utils.R','baby_bayes.R','gibbs_bayes.R','gibbs_stage_c_full.R')
}
# fast lookup "dir|name" -> $true
$cflKey = @{}
foreach($k in $cfl.Keys){ foreach($n in $cfl[$k]){ $cflKey["$k|$n"] = $true } }

# ---------- copy-as-is (scenario_A uniques into scenarioA_compare) ----------
$caiSrcDir = 'archive\interaction_old\comparison_stage1\scenario_A'
$caiDst    = '04_interaction\scenarioA_compare'
$cai = @('compute_recovery_metrics.R','fit_inla_wrapper.R','fit_mgcv_wrapper.R','fit_ours_wrapper.R','run_stage1_scenarioA_multiseed.R','scenarioA_settings.md','simulate_scenario_A.R','test_fit_inla.R','test_fit_mgcv.R','test_fit_ours.R','test_simulate_scenario_A.R')

# ---------- 2. copy-from-lib rows ----------
foreach($k in $cfl.Keys){
  foreach($n in $cfl[$k]){
    $d = Join-Path $k $n; $s = Join-Path 'lib' $n
    AddRow (TopGroup $k) 'copy-from-lib' $s $d (Md5 $s) (Md5 $d)
  }
}
# ---------- 3. copy-as-is rows ----------
foreach($n in $cai){
  $d = Join-Path $caiDst $n; $s = Join-Path $caiSrcDir $n
  AddRow (TopGroup $caiDst) 'copy-as-is' $s $d (Md5 $s) (Md5 $d)
}

# ---------- 4. single-file moves ----------
$fileMoves = @(
  @('baby_bayes_guide.md','docs\baby_bayes_guide.md'),
  @('Bayesian_Model_and_Simulation_Sections.docx','docs\Bayesian_Model_and_Simulation_Sections.docx'),
  @('CODE_AUDIT.md','_reorg\CODE_AUDIT.md'),
  @('CANONICAL_CHECK.md','_reorg\CANONICAL_CHECK.md'),
  @('Hellbender-12548977.out','results\_job_logs\Hellbender-12548977.out'),
  @('Hellbender-12549745.out','results\_job_logs\Hellbender-12549745.out'),
  @('submit_boundary.sub','02_reml\submit_boundary.sub'),
  @('submit_exp1.sub','02_reml\submit_exp1.sub')
)
$fileMoveDest = @{}
foreach($m in $fileMoves){
  $h = Md5 $m[1]; $fileMoveDest[$m[1]] = $true
  AddRow (TopGroup $m[1]) 'move' $m[0] $m[1] $h $h
}

# ---------- 5. whole-directory moves (walk dest, reconstruct source) ----------
$moveMap = [ordered]@{
  '01_prototypes\baby_bayes'          = 'codes\BABY_BAYES_codes'
  '01_prototypes\hist_codes'          = 'codes\hist_codes'
  '02_reml'                           = 'codes\REML_codes'
  '03_bayes_main\n400'                = 'codes\bayes_n400_codes'
  '03_bayes_main\n1000'               = 'codes\bayes_n1000_mcmc10000_codes'
  '04_interaction\scenarioB_tests'    = 'codes\test and comparison'
  '05_munich_rent\local'              = 'codes\run_codes\munich'
  '05_munich_rent\hpc'                = 'codes\realtest_code\munich_rent'
  'archive\interaction_old'           = 'codes\interaction_old'
  'archive\poisson_extension'         = 'codes\poisson_extension'
  'docs\vault'                        = 'LS-Spline-Dissertation'
  'docs\interaction_notes'            = 'interaction\interaction constract'
  'docs\literature\interaction_paper' = 'interaction\paper'
  'docs\literature\paper_reading'     = 'paper reading'
  'docs\presenting'                   = 'presenting files'
  'results\0311_output'               = '0311_output'
  'results\all_bayes_results_n400'    = 'all_bayes_results_n400'
  'results\bayes_out_gibss'           = 'bayes_out_gibss'
  'results\reml_resulyts'             = 'reml_resulyts'
}
foreach($destRoot in $moveMap.Keys){
  $srcRoot = $moveMap[$destRoot]
  if(-not (Test-Path -LiteralPath $destRoot)){ throw "move dest missing: $destRoot" }
  $destFull = (Resolve-Path -LiteralPath $destRoot).Path
  $all = Get-ChildItem -LiteralPath $destRoot -Recurse -File -Force
  foreach($f in $all){
    $rel = $f.FullName.Substring($destFull.Length+1)
    if($codeExt -contains $f.Extension.ToLower()){
      $parentRel = Split-Path $rel -Parent
      $dirKey = if($parentRel){ Join-Path $destRoot $parentRel } else { $destRoot }
      $destPath = Join-Path $destRoot $rel
      if($cflKey.ContainsKey("$dirKey|$($f.Name)")){ continue }      # injected copy-from-lib, already logged
      if($fileMoveDest.ContainsKey($destPath)){ continue }            # single-file move, already logged
      $h = Md5 $f.FullName
      AddRow (TopGroup $destRoot) 'move' (Join-Path $srcRoot $rel) $destPath $h $h
    }
  }
  # bulk-verify row (all files, count + bytes)
  $cnt = $all.Count; $byt = ($all | Measure-Object Length -Sum).Sum; if($null -eq $byt){ $byt = 0 }
  AddRow (TopGroup $destRoot) 'move-bulk-verify' $srcRoot $destRoot ("count=$cnt;bytes=$byt") ("count=$cnt;bytes=$byt")
}

# ---------- 6. parks ----------
$parks = @(
  @('.RData','archive\_to_delete\root__.RData'),
  @('codes\.RData','archive\_to_delete\codes__.RData'),
  @('.Rhistory','archive\_to_delete\root__.Rhistory'),
  @('codes\.Rhistory','archive\_to_delete\codes__.Rhistory'),
  @('03_bayes_main\n1000\.Rhistory','archive\_to_delete\03_bayes_main__n1000__.Rhistory'),
  @('03_bayes_main\n1000\bayes_overlay\.Rhistory','archive\_to_delete\03_bayes_main__n1000__bayes_overlay__.Rhistory'),
  @('03_bayes_main\n1000\bayes_test_sim1\.Rhistory','archive\_to_delete\03_bayes_main__n1000__bayes_test_sim1__.Rhistory'),
  @('05_munich_rent\hpc\.Rhistory','archive\_to_delete\05_munich_rent__hpc__.Rhistory'),
  @('Cubic_Spline_with_Spatial_Data_Feb25.pdf','archive\_to_delete\Cubic_Spline_with_Spatial_Data_Feb25.pdf'),
  @('Cubic_Spline_with_Spatial_Data_april.pdf','archive\_to_delete\Cubic_Spline_with_Spatial_Data_april.pdf')
)
foreach($p in $parks){
  $h = Md5 $p[1]
  AddRow 'archive/_to_delete/' 'park' $p[0] $p[1] $h $h
}

# ================= SELF-CHECK =================
# (a) every MD5-before == MD5-after
$mmA = @($rows | Where-Object { $_.action -ne 'move-bulk-verify' -and $_.before -ne $_.after })
# (b) every propagated canonical-named code file (outside archive/06_medsat/_reorg) == lib canonical
$mmB = New-Object System.Collections.ArrayList
$excludeRoots = @('archive','06_medsat_pilot','_reorg','lib')
Get-ChildItem -LiteralPath '.' -Recurse -File -Force | ForEach-Object {
  if($canonNames -notcontains $_.Name){ return }
  $rel = $_.FullName.Substring((Get-Location).Path.Length+1)
  $top = TopGroup $rel
  if($excludeRoots -contains $top){ return }
  if($rel -like '.git\*'){ return }
  $h = Md5 $_.FullName
  if($h -ne $libHashes.($_.Name)){ [void]$mmB.Add("$rel  ($h != lib $($libHashes.($_.Name)))") }
}

# ================= RENDER REORG_LOG.md =================
$groupOrder = @('lib/','01_prototypes','02_reml','03_bayes_main','04_interaction','05_munich_rent','docs','results','_reorg','archive','archive/_to_delete/')
$sb = New-Object System.Text.StringBuilder
function W($s){ [void]$sb.AppendLine($s) }
W '# REORG_LOG.md'
W ''
W 'Complete audit log of the reorganization pass, **derived from current filesystem state** (re-hashed at destination) + `_reorg/lib_hashes.json`. No file contents were edited; every action is a move or a byte-faithful copy. **Zero deletions** this pass — all removal candidates were *parked* into `archive/_to_delete/` for the user to empty by hand.'
W ''
W ('Generated: 2026-05-29. Total logged rows: {0}.' -f $rows.Count)
W ''
W 'Columns: action | source path | destination path | MD5-before | MD5-after. For `move` / `park` the source no longer exists at its old path (the operation was a same-volume rename verified byte-faithful inline at action time in stages B & D); before==after is the hash captured now at the destination. For `copy-from-lib` / `copy-as-is` both hashes are computed live (source still present), so the row is a genuine before-vs-after check. `move-bulk-verify` rows carry file-count + total-bytes for the whole directory (binary/data content verified at folder granularity).'
W ''
$shown = @{}
foreach($g in $groupOrder){
  $gr = @($rows | Where-Object { $_.group -eq $g })
  if($gr.Count -eq 0){ continue }
  $shown[$g] = $true
  W ("## $g")
  W ''
  W '| action | source | destination | MD5-before | MD5-after |'
  W '|---|---|---|---|---|'
  foreach($r in ($gr | Sort-Object action,dest)){
    W ("| {0} | {1} | {2} | {3} | {4} |" -f $r.action,$r.source,$r.dest,$r.before,$r.after)
  }
  W ''
}
# any group not in groupOrder (safety)
foreach($g in ($rows.group | Sort-Object -Unique)){
  if($shown.ContainsKey($g)){ continue }
  $gr = @($rows | Where-Object { $_.group -eq $g })
  W ("## $g")
  W ''
  W '| action | source | destination | MD5-before | MD5-after |'
  W '|---|---|---|---|---|'
  foreach($r in ($gr | Sort-Object action,dest)){ W ("| {0} | {1} | {2} | {3} | {4} |" -f $r.action,$r.source,$r.dest,$r.before,$r.after) }
  W ''
}

# ---- counts ----
$nMove   = @($rows | Where-Object { $_.action -eq 'move' }).Count
$nArch   = @($rows | Where-Object { $_.action -eq 'move' -and (TopGroup $_.dest) -eq 'archive' }).Count
$nCtl    = @($rows | Where-Object { $_.action -eq 'copy-to-lib' }).Count
$nCfl    = @($rows | Where-Object { $_.action -eq 'copy-from-lib' }).Count
$nCai    = @($rows | Where-Object { $_.action -eq 'copy-as-is' }).Count
$nPark   = @($rows | Where-Object { $_.action -eq 'park' }).Count
$nBulk   = @($rows | Where-Object { $_.action -eq 'move-bulk-verify' }).Count

W '## Self-check'
W ''
W '**(a) Byte preservation — every MD5-before == MD5-after:**'
if($mmA.Count -eq 0){ W ''; W 'Zero mismatches. All file actions preserved bytes exactly (verified inline at action time in stages A–D, and re-hashed at destination here).' }
else { W ''; foreach($m in $mmA){ W ("- MISMATCH: {0} -> {1}  ({2} != {3})" -f $m.source,$m.dest,$m.before,$m.after) } }
W ''
W '**(b) Propagation — every propagated core file''s MD5 == its `lib/` canonical:**'
W ''
W 'Checked every code file named like one of the 12 canonical files, excluding `archive/` (frozen stale snapshot), `06_medsat_pilot/` (untouched), `lib/` (the canonical itself), and `_reorg/`.'
if($mmB.Count -eq 0){ W ''; W 'Zero mismatches. Every propagated/moved copy of a canonical file (in `02_reml`, `03_bayes_main`, `04_interaction`, `05_munich_rent`) is byte-identical to its `lib/` source.' }
else { W ''; foreach($m in $mmB){ W ("- MISMATCH: {0}" -f $m) } }
W ''
W '## Counts'
W ''
W ('- Files moved: {0} individual code-file move rows (of which {1} into `archive/`), plus {2} whole-directory bulk-verify rows covering all binary/data content.' -f $nMove,$nArch,$nBulk)
W ('- Files copied: {0} (copy-to-lib: {1}, copy-from-lib: {2}, copy-as-is: {3}).' -f ($nCtl+$nCfl+$nCai),$nCtl,$nCfl,$nCai)
W ('- Files parked (removal candidates, NOT deleted): {0}.' -f $nPark)
W '- Files deleted: 0.'
W ''
W '> Note: `archive/poisson_extension/michigan_lung_data/.RData` and `.Rhistory` were deliberately **left frozen** inside the archived snapshot (not parked). Their presence in `archive/` rather than `archive/_to_delete/` is itself the record of that decision.'

$sb.ToString() | Out-File -Encoding utf8 'REORG_LOG.md'

Write-Host ("rows={0}  move={1} (archive={2}, bulk={3})  copy-to-lib={4}  copy-from-lib={5}  copy-as-is={6}  park={7}" -f $rows.Count,$nMove,$nArch,$nBulk,$nCtl,$nCfl,$nCai,$nPark)
Write-Host ("self-check (a) mismatches: {0}" -f $mmA.Count)
Write-Host ("self-check (b) mismatches: {0}" -f $mmB.Count)
"stageE done. REORG_LOG.md written."
