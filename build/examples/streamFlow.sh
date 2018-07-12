
psi="cos(2*pi*x)/(4*pi) + y**2"

# generate stream and vector data at higher resolution
python streamFlow.py --dims "101, 201" --psi "$psi" >& res.txt
mv streamFlow.vtk streamFlow101x201.vtk

# A
python streamFlow.py --xline "0.5 + t" --yline "0*t" --n 1 --psi "$psi" >& res.txt
exactA1=$(cat res.txt | grep integrated | awk '{print $5}')
errorA1=$(cat res.txt | grep integrated | awk '{print $8}')
mv streamLine.vtk streamLineA1.vtk

# B
python streamFlow.py --xline "0.8 + 0.6*cos(2*pi*t)" --yline "0.1 + 0.3*sin(2*pi*t)" --n 4 --psi "$psi" >& res.txt
exactB4=$(cat res.txt | grep integrated | awk '{print $5}')
errorB4=$(cat res.txt | grep integrated | awk '{print $8}')
python streamFlow.py --xline "0.8 + 0.6*cos(2*pi*t)" --yline "0.1 + 0.3*sin(2*pi*t)" --n 32 --psi "$psi" >& res.txt
exactB32=$(cat res.txt | grep integrated | awk '{print $5}')
errorB32=$(cat res.txt | grep integrated | awk '{print $8}')
mv streamLine.vtk streamLineB32.vtk

# C

python streamFlow.py --xline "0.2 + 1.6*sin(pi*t/2.)" --yline "-0.3 + 0.5*cos(pi*t/2.)" --n 1 --psi "$psi" >& res.txt
exactC1=$(cat res.txt | grep integrated | awk '{print $5}')
errorC1=$(cat res.txt | grep integrated | awk '{print $8}')

python streamFlow.py --xline "0.2 + 1.6*sin(pi*t/2.)" --yline "-0.3 + 0.5*cos(pi*t/2.)" --n 16 --psi "$psi" >& res.txt
exactC16=$(cat res.txt | grep integrated | awk '{print $5}')
errorC16=$(cat res.txt | grep integrated | awk '{print $8}')
mv streamLine.vtk streamLineC16.vtk

# D
python streamFlow.py --xline "0.25 + 1.6*sin(pi*t/2.)" --yline "-0.3 + 0.5*cos(pi*t/2.)" --n 1 --psi "$psi" >& res.txt
exactD1=$(cat res.txt | grep integrated | awk '{print $5}')
errorD1=$(cat res.txt | grep integrated | awk '{print $8}')

python streamFlow.py --xline "0.25 + 1.6*sin(pi*t/2.)" --yline "-0.3 + 0.5*cos(pi*t/2.)" --n 16 --psi "$psi" >& res.txt
exactD16=$(cat res.txt | grep integrated | awk '{print $5}')
errorD16=$(cat res.txt | grep integrated | awk '{print $8}')
mv streamLine.vtk streamLineD16.vtk

echo "A  1: exact=$exactA1  error=$errorA1"
echo "B  4: exact=$exactB4  error=$errorB4"
echo "B 32: exact=$exactB32 error=$errorB32"
echo "C  1: exact=$exactC1  error=$errorC1"
echo "C 16: exact=$exactC16 error=$errorC16"
echo "D  1: exact=$exactD1  error=$errorD1"
echo "D 16: exact=$exactD16 error=$errorD16"
