# Don't forget to check: https://tutorials.pumas.ai/html/DataWranglingInJulia/01-intro.html
# Additionally: https://tutorials.pumas.ai/html/DataWranglingInJulia/02-julia-syntax.html
# Also: https://tutorials.pumas.ai/html/DataWranglingInJulia/03-functions.html

##########################
#      Variable types    #
##########################

x = 5
typeof(x)

y = 3.14
typeof(y)

my_string = "Julia is Awesome"
typeof(my_string)

##########################
#        Containers      #
##########################

# Arrays - Vectors
typeof([1, 2, 3])
x_vec = [1, 2, 3]
x_vec_str = ["Hello", "Bye"]

# Arrays - Matrices
typeof([1 2 3])
x_matrix =[1 2 3; 4 5 6]

# Arrays are mutable
x_vec[2] = 42
x_matrix[1, 2] = 42
x_matrix[:, 1] = [55, 56]

# Tuples
x_tuple = (1,2,3)

# Tuple are immutable
x_tuple[1] = 42

# Tuples are not so convenient
user_options = ("Jose", 34, "Brazil")
# similar to R's lists: list("Jose", 34, "Brazil")

# NamedTuple
nt_options = (; age = 34, loc = "Brazil", name = "Jose")
# similar to R's named lists: list(age = 34, loc = "Brazil", name = "Jose")
nt_options.age

##########################
#         Functions      #
##########################

function fun_name(arg1, arg2)
    some_computation = arg1 ∘ arg2 # or some other crazy some_computation
end

# Compact assignment form
round_number(x) = round(x)

# it works with any type
round_number(3.14)
round_number(5)

# Specialized methods and a quick intro into multiple dispatch
round_number(x::Float64) = round(x)   # this is not a type "annotation", this is way more strict
round_number(x::Int64) = x

@which round_number(3.14)
@which round_number(5)

# Take a look at one example in the documentation for `simobs`: https://docs.pumas.ai/stable/pumas_docstrings/#Pumas.simobs

##########################
#  Anonymous Functions   #
##########################

typeof(x -> x^2)  # simple reference

map(x -> x^2, 3)                  # single value
map(x -> x^2, 1:10)               # collection of value
map((x, y) -> x^2 + y, 2, 1)      # more values as Tuples
map((x, y) -> x^2 + y, 1:3, 3:5)  # more values as Tuples applied to collections

# I can also make any anymous function a function (no, it is not that obvious)
(x -> x^2)(3)
my_anon_function = x -> x^2
my_anon_function(3)

# Example in Pumas population construction

map(id -> Subject(; id=id), )

##########################
#  Array comprehensions  #
##########################

# Similar to Python's list comprehensions

[x^2 for x in 1:10]

# multiple inputs
[x * y for x in 1:10 for y in 1:2]

# conditionals
[x^2 for x in 1:10 if isodd(x)]

# using a desired Type
Float64[x^2 for x in 1:10 if isodd(x)]
Complex[x^2 for x in 1:10 if iseven(x)]

##########################
#  Ternary operator `?`  #
##########################

# Syntax: condition ? value_if_true : value_if_false

var_1 = iseven(2) ? "yes" : "no"
var_2 = isodd(2) ? "yes" : "no"
