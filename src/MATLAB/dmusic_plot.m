function dmusic_plot(pfs, fs, ts)

		figure;
		imagesc(ts, fs, pfs)
		set(gca, 'Ydir', 'normal')
		xlabel 'Time [ms]'
		ylabel 'Frequency [kHz]'
		colorbar

end