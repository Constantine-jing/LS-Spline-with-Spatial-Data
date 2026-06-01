$ErrorActionPreference = 'Stop'
function Md5($p){ (Get-FileHash -Algorithm MD5 -LiteralPath $p).Hash }
$codeExt = @('.r','.cpp','.h','.sub','.md','.sh')
$root = (Get-Location).Path

$log = New-Object System.Collections.ArrayList
@(Get-Content '_reorg\reorg_actions.json' -Raw | ConvertFrom-Json) | ForEach-Object { [void]$log.Add($_) }

function Move-Tree($src,$dst,$group){
  if(-not (Test-Path -LiteralPath $src)){ Write-Host "SKIP (missing): $src"; return }
  $srcFull = (Resolve-Path -LiteralPath $src).Path
  # pre: code files keyed by relpath
  $preCode = @{}
  Get-ChildItem -LiteralPath $src -Recurse -File -Force | ForEach-Object {
    if($codeExt -contains $_.Extension.ToLower()){
      $rel = $_.FullName.Substring($srcFull.Length+1)
      $preCode[$rel] = (Md5 $_.FullName)
    }
  }
  # pre: bulk stats (all files)
  $preAll = Get-ChildItem -LiteralPath $src -Recurse -File -Force
  $preCount = $preAll.Count
  $preBytes = ($preAll | Measure-Object Length -Sum).Sum
  if($null -eq $preBytes){ $preBytes = 0 }
  # ensure dest, move children in
  if(-not (Test-Path -LiteralPath $dst)){ New-Item -ItemType Directory -Force -Path $dst | Out-Null }
  Get-ChildItem -LiteralPath $src -Force | Move-Item -Destination $dst
  $dstFull = (Resolve-Path -LiteralPath $dst).Path
  # post: per-code-file verify + log
  foreach($rel in $preCode.Keys){
    $dp = Join-Path $dstFull $rel
    $a = Md5 $dp
    if($a -ne $preCode[$rel]){ throw "MOVE HASH MISMATCH: $rel ($src -> $dst)" }
    [void]$log.Add([ordered]@{ group=$group; action='move'; source=(Join-Path $src $rel); dest=$dp; before=$preCode[$rel]; after=$a })
  }
  # post: bulk verify
  $postAll = Get-ChildItem -LiteralPath $dst -Recurse -File -Force
  $postCount = $postAll.Count
  $postBytes = ($postAll | Measure-Object Length -Sum).Sum
  if($null -eq $postBytes){ $postBytes = 0 }
  if($postCount -ne $preCount -or $postBytes -ne $preBytes){ throw "BULK MISMATCH $src -> $dst : count $preCount/$postCount bytes $preBytes/$postBytes" }
  [void]$log.Add([ordered]@{ group=$group; action='move-bulk-verify'; source="$src ($preCount files, $preBytes bytes)"; dest="$dst ($postCount files, $postBytes bytes)"; before="count=$preCount;bytes=$preBytes"; after="count=$postCount;bytes=$postBytes" })
  # remove emptied source
  if((Get-ChildItem -LiteralPath $src -Force | Measure-Object).Count -eq 0){ Remove-Item -LiteralPath $src -Force }
  Write-Host ("MOVED {0} -> {1}  [code={2} bulk={3} files]" -f $src,$dst,$preCode.Count,$preCount)
}

Move-Tree 'codes\BABY_BAYES_codes'              '01_prototypes\baby_bayes'      '01_prototypes/baby_bayes/'
Move-Tree 'codes\hist_codes'                    '01_prototypes\hist_codes'      '01_prototypes/hist_codes/'
Move-Tree 'codes\REML_codes'                    '02_reml'                       '02_reml/'
Move-Tree 'codes\bayes_n400_codes'              '03_bayes_main\n400'            '03_bayes_main/n400/'
Move-Tree 'codes\bayes_n1000_mcmc10000_codes'   '03_bayes_main\n1000'           '03_bayes_main/n1000/'
Move-Tree 'codes\test and comparison'           '04_interaction\scenarioB_tests' '04_interaction/scenarioB_tests/'
Move-Tree 'codes\run_codes\munich'              '05_munich_rent\local'          '05_munich_rent/local/'
Move-Tree 'codes\realtest_code\munich_rent'     '05_munich_rent\hpc'            '05_munich_rent/hpc/'
Move-Tree 'codes\interaction_old'               'archive\interaction_old'       'archive/interaction_old/'
Move-Tree 'codes\poisson_extension'             'archive\poisson_extension'     'archive/poisson_extension/'

# remove now-empty intermediate shells
foreach($e in @('codes\run_codes','codes\realtest_code')){
  if((Test-Path -LiteralPath $e) -and ((Get-ChildItem -LiteralPath $e -Force | Measure-Object).Count -eq 0)){ Remove-Item -LiteralPath $e -Force }
}

($log | ConvertTo-Json -Depth 5) | Out-File -Encoding utf8 '_reorg\reorg_actions.json'
"stageB done. total log rows: $($log.Count)"
