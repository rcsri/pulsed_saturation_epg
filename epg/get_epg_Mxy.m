function Mxy = get_epg_Mxy(FZ, flip_phase)

Mxy = abs( FZ(1,1) * exp(-1i * flip_phase)); % transverse (0th state), Phase-Demodulated signal.
