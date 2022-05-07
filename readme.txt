
This is a routine for analysis of Nortek Vectrino ADV data as used in the fluvial hydraulics lab at UBC Vancouver.

The ADV is mounted to a ZaberMotion linear stage. 
The stage is programmed by "stage_mover.py" to iterate through a sequence of positions as the ADV collects data.
This script will need to be customized for each application.

The flow data are saved from the Nortek Vectrino software with a base filename and a timestamp.
The timestamp in the filename is important since it allows correlation between the stage location and velocity timeseries.

Functionality to obtain flow statistics are contained in "adv_analysis.py". 
Data are despiked prior to analysis using the method of Nikora and Goring 1998.

Example data is contained in the "example-data" directory.

Example usage of the routine is provided in "example_usage.py".
I.e., try "python example_usage.py", fix the missing dependencies, then try again
The analysis output will be generated in the "example_data" directory.

Feel free to direct questions to kpierce@alumni.ubc.ca
