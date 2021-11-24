
For ($i=1; $i -le 9; $i++) {
    $file = 'C:\Users\amcga\GitHub\episim\simulation_results\full_run_' + $i + '.pickle'

    if (-not(Test-Path -Path $file -PathType Leaf)) {
        Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", "C:\Users\amcga\GitHub\episim\sims_$i.ipynb"
    }
}

<#


Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", 'C:\Users\amcga\GitHub\episim\sims_2.ipynb'
Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", 'C:\Users\amcga\GitHub\episim\sims_3.ipynb'
Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", 'C:\Users\amcga\GitHub\episim\sims_4.ipynb'
Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", 'C:\Users\amcga\GitHub\episim\sims_5.ipynb'
Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", 'C:\Users\amcga\GitHub\episim\sims_6.ipynb'
Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", 'C:\Users\amcga\GitHub\episim\sims_7.ipynb'
Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", 'C:\Users\amcga\GitHub\episim\sims_8.ipynb'
Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute","--ExecutePreprocessor.timeout=180000", 'C:\Users\amcga\GitHub\episim\sims_9.ipynb'

$jobs = @();

$j = Start-Job -ScriptBlock {
    jupyter-nbconvert --inplace --execute '2 exec sims 1.ipynb'
}
$jobs.Add($j)

Start-Process -NoNewWindow "C:\Users\amcga\envs\citation-deaths\Scripts\jupyter-nbconvert.exe" -ArgumentList "--inplace", "--execute", 'C:\Users\amcga\GitHub\episim\2/ exec/ sims/ 1.ipynb'

$j = Start-Job -ScriptBlock {
    jupyter-nbconvert --inplace --execute '2 exec sims 3.ipynb';
    jupyter-nbconvert --inplace --execute '2 exec sims 4.ipynb';
}
$jobs.Add($j)

$j = Start-Job -ScriptBlock {
    jupyter-nbconvert --inplace --execute '2 exec sims 5.ipynb';
    jupyter-nbconvert --inplace --execute '2 exec sims 6.ipynb';
}
$jobs.Add($j)

$j = Start-Job -ScriptBlock {
    jupyter-nbconvert --inplace --execute '2 exec sims 7.ipynb';
    jupyter-nbconvert --inplace --execute '2 exec sims 8.ipynb';
}
$jobs.Add($j)

for ($i=0; $i -lt $jobs.Count; $i++) {
    wait-job $jobs[$i];
}
#>
<#
For ($i=1; $i -le 10; $i++) {

    $jstart = {
        param($j)
        Write-Host "HELLO";
        jupyter-nbconvert --inplace --execute '2 exec sims $j.ipynb';
    }

    Start-Job $jstart -Arg $i
    Write-Host jupyter-nbconvert --inplace --execute '2 exec sims $i.ipynb';

    $jobs.Add( $job );
}

for ($i=0; $i -lt $jobs.Count; $i++) {
    wait-job $jobs[$i];
    Write-Host 
}
#>