$ErrorActionPreference = 'Stop'
function Md5($p){ (Get-FileHash -Algorithm MD5 -LiteralPath $p).Hash }

# Move a whole directory by RENAME to an explicit final path (parent must exist, final must NOT).
# Verified by bulk file-count + total-bytes (same-volume rename cannot alter bytes).
function MoveDir($src,$finalDst){
  if(-not (Test-Path -LiteralPath $src)){ throw "MISSING src: $src" }
  if(Test-Path -LiteralPath $finalDst){ throw "UNEXPECTED existing dest: $finalDst" }
  $pre = Get-ChildItem -LiteralPath $src -Recurse -File -Force
  $pc = $pre.Count; $pb = ($pre | Measure-Object Length -Sum).Sum; if($null -eq $pb){ $pb = 0 }
  Move-Item -LiteralPath $src -Destination $finalDst
  $post = Get-ChildItem -LiteralPath $finalDst -Recurse -File -Force
  $qc = $post.Count; $qb = ($post | Measure-Object Length -Sum).Sum; if($null -eq $qb){ $qb = 0 }
  if($qc -ne $pc -or $qb -ne $pb){ throw "DIR BULK MISMATCH $src -> $finalDst : count $pc/$qc bytes $pb/$qb" }
  Write-Host ("DIR  {0} -> {1}  [{2} files, {3} bytes]" -f $src,$finalDst,$qc,$qb)
}

# Move a single file to an explicit final path; verify MD5 before == after.
function MoveFile($src,$finalDst){
  if(-not (Test-Path -LiteralPath $src)){ throw "MISSING src: $src" }
  if(Test-Path -LiteralPath $finalDst){ throw "UNEXPECTED existing dest: $finalDst" }
  $b = Md5 $src
  Move-Item -LiteralPath $src -Destination $finalDst
  $a = Md5 $finalDst
  if($a -ne $b){ throw "FILE HASH MISMATCH $src -> $finalDst" }
  Write-Host ("FILE {0} -> {1}  [{2}]" -f $src,$finalDst,$a)
}

# Park a removal candidate into archive\_to_delete\ under a collision-safe name; verify MD5.
function Park($src,$newName){
  $dst = Join-Path 'archive\_to_delete' $newName
  if(-not (Test-Path -LiteralPath $src)){ throw "MISSING src: $src" }
  if(Test-Path -LiteralPath $dst){ throw "UNEXPECTED existing dest: $dst" }
  $b = Md5 $src
  Move-Item -LiteralPath $src -Destination $dst
  $a = Md5 $dst
  if($a -ne $b){ throw "PARK HASH MISMATCH $src -> $dst" }
  Write-Host ("PARK {0} -> {1}  [{2}]" -f $src,$dst,$a)
}

# ---- DOCS ----
MoveDir 'LS-Spline-Dissertation'        'docs\vault'
MoveDir 'interaction\interaction constract' 'docs\interaction_notes'
MoveDir 'interaction\paper'             'docs\literature\interaction_paper'
MoveDir 'paper reading'                 'docs\literature\paper_reading'
# docs\presenting pre-exists empty: assert empty, remove it, then rename source in.
if((Get-ChildItem -LiteralPath 'docs\presenting' -Force | Measure-Object).Count -ne 0){ throw "docs\presenting not empty" }
Remove-Item -LiteralPath 'docs\presenting' -Force
MoveDir 'presenting files'              'docs\presenting'
MoveFile 'baby_bayes_guide.md'          'docs\baby_bayes_guide.md'
MoveFile 'Bayesian_Model_and_Simulation_Sections.docx' 'docs\Bayesian_Model_and_Simulation_Sections.docx'
MoveFile 'CODE_AUDIT.md'                '_reorg\CODE_AUDIT.md'
MoveFile 'CANONICAL_CHECK.md'           '_reorg\CANONICAL_CHECK.md'

# ---- RESULTS ----
MoveDir '0311_output'                   'results\0311_output'
MoveDir 'all_bayes_results_n400'        'results\all_bayes_results_n400'
MoveDir 'bayes_out_gibss'               'results\bayes_out_gibss'
MoveDir 'reml_resulyts'                 'results\reml_resulyts'
MoveFile 'Hellbender-12548977.out'      'results\_job_logs\Hellbender-12548977.out'
MoveFile 'Hellbender-12549745.out'      'results\_job_logs\Hellbender-12549745.out'

# ---- PARK (no deletions; user empties archive\_to_delete\ by hand) ----
Park '.RData'                                              'root__.RData'
Park 'codes\.RData'                                        'codes__.RData'
Park '.Rhistory'                                           'root__.Rhistory'
Park 'codes\.Rhistory'                                     'codes__.Rhistory'
Park '03_bayes_main\n1000\.Rhistory'                       '03_bayes_main__n1000__.Rhistory'
Park '03_bayes_main\n1000\bayes_overlay\.Rhistory'         '03_bayes_main__n1000__bayes_overlay__.Rhistory'
Park '03_bayes_main\n1000\bayes_test_sim1\.Rhistory'       '03_bayes_main__n1000__bayes_test_sim1__.Rhistory'
Park '05_munich_rent\hpc\.Rhistory'                        '05_munich_rent__hpc__.Rhistory'
Park 'Cubic_Spline_with_Spatial_Data_Feb25.pdf'           'Cubic_Spline_with_Spatial_Data_Feb25.pdf'
Park 'Cubic_Spline_with_Spatial_Data_april.pdf'           'Cubic_Spline_with_Spatial_Data_april.pdf'

# ---- CLEANUP EMPTY SHELLS ----
foreach($e in @('codes','interaction')){
  if(Test-Path -LiteralPath $e){
    $n = (Get-ChildItem -LiteralPath $e -Force | Measure-Object).Count
    if($n -eq 0){ Remove-Item -LiteralPath $e -Force; Write-Host "RMDIR $e (empty)" }
    else { throw "NOT EMPTY, not removed: $e ($n items)" }
  }
}

"stageD done."
