#!/bin/bash
#Youtube
grec -g youtube/com-youtube.ungraph.bin -P -o results_youtube.txt -N 10
grec -g youtube/com-youtube.ungraph.bin -D -o results_youtube.txt -N 10
grec -g youtube/com-youtube.ungraph.bin -E -o results_youtube.txt -N 10
#Livejournal
grec -g livejournal/com-lj.ungraph.bin -P -o results_livejournal.txt -N 10
grec -g livejournal/com-lj.ungraph.bin -D -o results_livejournal.txt -N 10
grec -g livejournal/com-lj.ungraph.bin -E -o results_livejournal.txt -N 10
#Orkut
#grec -g orkut/com-orkut.ungraph.bin -P -o results_orkut.txt -N 10
grec -g orkut/com-orkut.ungraph.bin -D -o results_orkut.txt -N 10
grec -g orkut/com-orkut.ungraph.bin -E -o results_orkut.txt -N 10
#Friendster
grec -g friendster/com-friendster.ungraph.bin -E -o results_friendster.txt -N 10
