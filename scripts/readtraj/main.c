/* Quick program for converting my measure files (in .txt)
* to Kari's format (binary with a header).
* Goal: use Kari's FSH program to reweight my measure files,
* so here I attempt to reproduce the I/O in Kari's main/measure.c
*/

#include "stuff.h"

int main(int argc, char *argv[]) {

	if (argc != 3) {
		printf("Usage: ./<program name> <trajectory file> <jackknife blocks>\n");
		die(0);
	}

	int jack_blocks = atoi(argv[2]);

	int sample = 1; // how often to 'sample' trajectories
	realtime(argv[1], sample, jack_blocks);

	return 0;
}
