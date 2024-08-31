grep {SP.* j01_setup_sector_mappings.log | sed 's/}/;/g' | sed 's/{//g' | sed 's/==/=/g' | sed 's/,SP/;\nSP/g' | sed 's/SP/id SP/g'
