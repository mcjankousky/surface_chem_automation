# /usr/bin/env python

import sys, os
import pandas as pd

calc_lst = []
E_lst = []
for root, dirs, files in os.walk("."):
	for name in dirs:
		paff = os.path.join(root, name)
		if os.path.exists(paff + '/OUTCAR'):
			calc_lst.append(paff)
			conv = False
			with open(paff + '/OUTCAR') as outcar:
				if 'Elaps' in outcar.read():
					conv = True
			#conv_lst.append(conv)	
			if conv == True:
				if os.path.exists(paff + '/OSZICAR'):
					with open(paff + '/OSZICAR') as oszicar:
						oszicar = oszicar.readlines()
						if len(oszicar) > 0:
							last_line = oszicar[-1]
							last_line = last_line.split()
							last_line = [x for x in last_line if x != '']
							Energy = last_line[2]
						else:
							Energy = 0
				else:
					Energy = 0
			else:
				Energy = 0
			E_lst.append(Energy)
				
E_dict = {}
E_dict['calculation'] = site_lst
E_dict['E'] = E_lst

E_df = pd.DataFrame.from_dict(E_dict)

E_df = E_df[['calculation','E']]

E_df.to_csv(sys.argv[1] + '.csv')
