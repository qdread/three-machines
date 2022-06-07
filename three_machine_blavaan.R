library(blavaan)
library(data.table)

# Set file path
fp <- ifelse(Sys.info()['sysname'] == 'Windows', 'soybean_machine_comparison/project', '/project/qdr/ars_misc_data')