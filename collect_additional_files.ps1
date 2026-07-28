# collect_additional_files.ps1
# Gathers every file cited as an Additional file into one folder, ready for
# upload to BMC Bioinformatics. Run from the project root:
#     .\collect_additional_files.ps1

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$dest = Join-Path $root "additional_files"
New-Item -ItemType Directory -Force -Path $dest | Out-Null

Write-Host "Collecting into $dest" -ForegroundColor Cyan
Write-Host ""

# Each entry: the Additional file number, what it is, and the patterns to look
# for anywhere under the project root.
$wanted = @(
    @{
        Num      = 1
        Title    = "Candidate spreadsheet"
        Patterns = @("PlantDefenseESM_Supplementary_Candidates*.xlsx")
        RenameTo = "PlantDefenseESM_Supplementary_Candidates.xlsx"
    },
    @{
        Num      = 2
        Title    = "Validation keyword partition (39 defense-specific + 8 broad)"
        Patterns = @("*keyword*partition*.csv", "*keyword*subset*.csv",
                     "*keyword*partition*.xlsx", "*validation_keywords*.csv")
        RenameTo = $null
    },
    @{
        Num      = 3
        Title    = "Threshold-sensitivity table"
        Patterns = @("*threshold*sensitivity*.csv", "*threshold*sensitivity*.xlsx",
                     "*sensitivity*threshold*.csv")
        RenameTo = $null
    },
    @{
        Num      = 4
        Title    = "Fully labelled top-50 heatmaps"
        Patterns = @("figS_category_heatmaps_full_*.png")
        RenameTo = $null
    }
)

$missing = @()

foreach ($item in $wanted) {
    $found = @()
    foreach ($pat in $item.Patterns) {
        $found += Get-ChildItem -Path $root -Recurse -File -Filter $pat `
                    -ErrorAction SilentlyContinue |
                  Where-Object { $_.DirectoryName -ne $dest }
    }
    $found = $found | Sort-Object FullName -Unique

    if ($found.Count -eq 0) {
        Write-Host ("Additional file {0}  [MISSING]  {1}" -f $item.Num, $item.Title) `
                   -ForegroundColor Yellow
        $missing += $item
        continue
    }

    foreach ($f in $found) {
        $name = if ($item.RenameTo -and $found.Count -eq 1) { $item.RenameTo } else { $f.Name }
        Copy-Item $f.FullName (Join-Path $dest $name) -Force
        $mb = [math]::Round($f.Length / 1MB, 2)
        Write-Host ("Additional file {0}  {1,-45} {2} MB" -f $item.Num, $name, $mb) `
                   -ForegroundColor Green
        if ($mb -gt 20) {
            Write-Host ("   WARNING: over the 20 MB per-file limit") -ForegroundColor Red
        }
    }
}

Write-Host ""
if ($missing.Count -gt 0) {
    Write-Host "Not found - these are cited in the manuscript but have no file:" `
               -ForegroundColor Yellow
    foreach ($m in $missing) {
        Write-Host ("   Additional file {0}: {1}" -f $m.Num, $m.Title)
    }
    Write-Host ""
}

Write-Host "Contents of $dest :" -ForegroundColor Cyan
Get-ChildItem $dest | Select-Object Name, @{N='MB';E={[math]::Round($_.Length/1MB,2)}}
