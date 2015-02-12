#!/usr/bin/env python
with open('berea_xsection_bot.dat', 'r') as fbot:
    bot_data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), fbot)
    v_solid = sum(map(lambda x: sum(x), bot_data))
    v_total = sum(map(lambda x: len(x), bot_data))
    print v_solid, v_total, float(v_solid)/v_total

with open('berea_xsection_top.dat', 'r') as ftop:
    top_data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), ftop)
    v_solid = sum(map(lambda x: sum(x), top_data))
    v_total = sum(map(lambda x: len(x), top_data))
    print v_solid, v_total, float(v_solid)/v_total
