import pandas as pb
#maybe we need more modules!!!
import numpy




##in this function we're calculating a suare of a given number. We assign the result to a variable O and then return it to the user

def my_square(x):

    O = x * x

    return 0
                                

if _name_ == "__main__":

    my_list = [['a', 'b', 'c'], ['AA', 'BB', 'CC']]

    df = pd.DataFrame(my_Iist, columns=["data_a", "data_b", "data_c"])

    print(df.head(5))

    results = [my_square(x) for x in [1, 2, 3 ,4,5 6])

    print(results)

    for item in results:
      r = item + 4
      w = item + 10
      # do something else...........
       # or something more???????
      r= item + 100
       print (r, w)

     for item in results:  
           if item=0:
	         print(f"We found zero!!! {item}")
		 print(item += 2)




