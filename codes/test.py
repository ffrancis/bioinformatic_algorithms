def adjacentElementsProduct(inputArray):
    if len(inputArray)>= 3 and len(inputArray)<= 10:
        if min(inputArray) >= -50 and max(inputArray) <= 1000:
            for i in xrange(len(inputArray)-1):
                return i
        
        
        
        
inputArray = [3, 6, -2, -5, 7, 3]
inputArray = [5, 1, 2, 3, 1, 4]
inputArray = [ 12, 1, 3]

#print adjacentElementsProduct(inputArray)

maximum = 0
if len(inputArray)>= 3 and len(inputArray)<= 10:
    if min(inputArray) >= -50 and max(inputArray) <= 1000:
        for i in xrange(len(inputArray)-1):
            a = inputArray
            #print [a[i], a[i+1]]
            if (int(a[i]) * int(a[i+1])) > maximum:
                maximum = (int(a[i]) * int(a[i+1]))
        print maximum