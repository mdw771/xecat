pl_noscat, = plt.plot(t, i_noscat_b, label='Unscattered', color='salmon')
pl_1el, = plt.plot(t, i_1el_b, label='Single elastically scattered', color='black')
pl_elpl, = plt.plot(t, i_elpl_b, label='Plural elastically scattered', color='darkcyan')
pl_inel, = plt.plot(t, i_inel_b, label='Inelastically scattered', color='magenta')
pl_out, = plt.plot(t, i_out_b, label='Scattered out', color='orange')
plt.semilogy()