Model ran from 2009-06-01 00:00:00 to 2009-09-01 00:00:00

Output from calibrated 1D model run including:
\output\temp_diff01.csv	-	time (rows) over depth (columns) matrix including temperature data after 1st modeling step (diffusion)
\output\temp_mix02.csv	-	time (rows) over depth (columns) matrix including temperature data after 2nd modeling step (mixing)
\output\temp_conv03.csv	-	time (rows) over depth (columns) matrix including temperature data after 3rd modeling step (convection)
\output\temp_total04.csv	-	time (rows) over depth (columns) matrix including temperature data after final (4th) modeling step (total)
\output\diff.csv	-	time (rows) over depth (columns) matrix including diffusion data before 1st modeling step (diffusion)
\output\buoyancy.csv	-	time (rows) over depth (columns) matrix including bouyancy data before 1st modeling step (diffusion)
\output\meteorology_input.csv	-	time (rows) over depth (columns) matrix including meteorology input data during all modeling steps

Output from calibrated 1D models GLM and Simstrat (both k-Epsilon):
\output\diff_gotm.csv	-	time (rows) over depth (columns) matrix including diffusivity data from GOTM model (diffusion)
\output\diff_simstrat.csv	-	time (rows) over depth (columns) matrix including diffusivity data from Simstrat model (diffusion)

Field data from under-water monitoring chain:
\output\observed_temp.csv	-	time (rows) over depth (columns) matrix of measured temperature data

Pseudo scheme:
temp_total04[t-dt] + meteorology_input[t] + buoyancy[t] --> diff[t] --> temp_diff01[t] --> temp_mix02[t] --> temp_conv03[t] --> temp_total04[t]
