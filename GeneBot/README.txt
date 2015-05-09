Run the quick-start shell(quick_start.sh) with recommended configuration
If you want to configure the bot manually, type the following command  for more configuration details:
python3 source/geneBot.py -h

After the analysis table is done, you can run the filter shell(filter.sh) with recommended configuration to extract proper guides.
Same as Above, If you want to configure the filter manually, type the following command  for more configuration details:
python3 source/tableFilter.py -h

Most errors occurs during waiting the analyzing server, so there is a recovery shell(recover.sh) that reads the backup file(analysis job backup.txt) and recover the lost analysis from the errors. Notice, if other reason causes the error, this method is not suitable. You should collect failed search words and run the Gene Bot again.
As usual, you still can configure the recover program. Type the following command for more configuration details:
python3 source/analysisFetcher.py -h 