$ErrorActionPreference = 'Stop'
function Md5($p){ (Get-FileHash -Algorithm MD5 -LiteralPath $p).Hash }

$munich = 'codes\run_codes\munich'
$tc = 'codes\test and comparison'
$srcMap = [ordered]@{
  'baby_bayes.R'                  = $munich
  'fit_spatial_reml.R'            = $munich
  'gibbs_bayes.R'                 = $munich
  'gibbs_interaction.R'           = $munich
  'gibbs_stage_c_full.R'          = $munich
  'ls_basis.R'                    = $munich
  'ls_interaction.R'              = $munich
  'ls_interaction_core.cpp'       = $munich
  'marginal_utils.R'              = $munich
  'spatial_utils.R'               = $munich
  'canonical_schema.R'            = $tc
  'fit_ours_interaction_wrapper.R'= $tc
}

$log = New-Object System.Collections.ArrayList
$libHashes = [ordered]@{}
foreach($f in $srcMap.Keys){
  $s = Join-Path $srcMap[$f] $f
  $d = Join-Path 'lib' $f
  $b = Md5 $s
  Copy-Item -LiteralPath $s -Destination $d
  $a = Md5 $d
  if($a -ne $b){ throw "lib copy mismatch for $f" }
  $libHashes[$f] = $a
  [void]$log.Add([ordered]@{ group='lib/'; action='copy-to-lib'; source=$s; dest=$d; before=$b; after=$a })
}
($log | ConvertTo-Json -Depth 5) | Out-File -Encoding utf8 '_reorg\reorg_actions.json'
($libHashes | ConvertTo-Json) | Out-File -Encoding utf8 '_reorg\lib_hashes.json'
"lib populated: $($libHashes.Count) files"
$libHashes.GetEnumerator() | ForEach-Object { "{0}  {1}" -f $_.Value, $_.Key }
