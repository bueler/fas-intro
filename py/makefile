# testing

runfas1_1:
	-@./testit.sh fas1.py "-monitor -monitorupdate -cyclemax 2 -coarse 20 -levels 2" 1 1

runfas1_2:
	-@./testit.sh fas1.py "-monitor -cyclemax 5 -rtol 0.0 -down 20 -ngsonly -mms" 1 2

runfas1_3:
	-@./testit.sh fas1.py "-monitor -cyclemax 3 -K 4 -down 2 -up 2 -mms -coarse 5" 1 3

runfas1_4:
	-@./testit.sh fas1.py "-K 5 -mms -fcycle -cyclemax 1" 1 4

runfas1_5:
	-@./testit.sh fas1.py "-monitor -K 4 -mms -fcycle -cyclemax 2 -R inj" 1 5

runfas1_6:
	-@./testit.sh fas1.py "-monitor -K 4 -mms -fcycle -cyclemax 2 -R inj -unroll" 1 6

test_fas1: runfas1_1 runfas1_2 runfas1_3 runfas1_4 runfas1_5 runfas1_6

test: test_fas1

.PHONY: clean runfas1_1 runfas1_2 runfas1_3 runfas1_4 runfas1_5 runfas1_6 test_fas1 test

clean:
	@rm -f maketmp tmp difftmp
	@rm -rf __pycache__
