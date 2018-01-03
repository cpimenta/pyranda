

# Copy and paste baselines here


baselines = """
euler-50 -- 9.65612670412e-09
euler-100 -- 7.6111541815e-11
euler-200 -- 5.93285138337e-13
euler-300 -- 3.47125615129e-14
euler-400 -- 4.63496183249e-15
"""





# Make dictionary
dbase = {}
for bb in baselines.split('\n'):
    if bb:
        name = bb.split('--')[0].strip()
        diff = bb.split('--')[1].strip()
        dbase[name] = diff
