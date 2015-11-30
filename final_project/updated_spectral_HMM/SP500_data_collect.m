% pulls in returns of SP500 close prices
[SP500, SP500_headers] = Quandl.get('YAHOO/INDEX_GSPC','start_date','2005-11-23','end_date','2015-11-23','transformation','rdiff');
data = SP500.Close.data;