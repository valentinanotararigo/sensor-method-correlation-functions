# sensor-method-correlation-functions

This code computes the power spectrum and the frequency- and time-resolved correlation functions of a two-level emitter in resonance with a cavity mode and coupled to it, using the sensor method. The system is described by a Linblad master equation, which involves incoherent pumping of the emitter and loss from the cavity and the atom. The analysis is carried out in the steady state condition.

The following figure represents the power spectra for three different values of the damping of the mode:

![power_spectra](power_spectrum/func_jk/ps_JC.pdf)

The figure depicts the sfrequency-filtered second order correlation function at zero time delay, again for three different values of the damping of the mode:

![g2(0)](second_order_g/frequency_resolved/g2_JC.pdf)

The following figure reports the time-resolved second order correlation function, where positive times correspond to detection of frequencies as reported in the legend, negative times to the opposite order:

![g2(tau)](second_order_g/time_resolved/g2t_JC.pdf)
