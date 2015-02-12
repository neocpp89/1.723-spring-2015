#!/usr/bin/env python
def calculate_porosity(data):
    v_solid = sum(map(lambda x: sum(x), data))
    v_total = sum(map(lambda x: len(x), data))
    return (v_solid, v_total, float(v_solid)/v_total)

if __name__ == "__main__":
    with open('berea_xsection_bot.dat', 'r') as fbot:
        bot_data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), fbot)
        print calculate_porosity(bot_data) 

    with open('berea_xsection_top.dat', 'r') as ftop:
        top_data = map(lambda x: map(lambda y: int(y), x.rstrip().split(' ')), ftop)
        print calculate_porosity(top_data) 
