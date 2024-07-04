release: clean compileRelease run plot
debug: clean compileDebug run plot
python: clean runPython plot

plot:
	./plotRoutes.gmt
	./plotConvergence.gmt

run:
	time ./tsp

compileRelease:
	sed -i 's/\t/    /g' tsp.F
	gfortran -O3 -ffree-form -std=f2018 -o tsp tsp.F
	chmod +x tsp

compileDebug:
	sed -i 's/\t/    /g' tsp.F
	gfortran -Ddebug -g -fcheck=bounds -ffree-form -std=f2018 -o tsp tsp.F
	chmod +x tsp

runPython:
	sed -i 's/    /\t/g' tsp.py
	time python3 tsp.py


clean:
	rm -f tsp
	rm -f *.dat *.pdf *.jpg
	rm -f gmt.conf
	rm -f *.mod
