
# Simulation distribution sum_{j=1}^n U_j(U_j+1)...(U_j+k-2)

B=10^5
m=10
n=50
vec = rwilcox(nn=B, m=m, n=n)
vec2 = vec^2
plot(table(vec2))

vec3 = vec^3
plot(table(vec3))

vec4 = vec^4
plot(table(vec4))

vec5 = vec^5
plot(table(vec5))

vec6 = vec^6
plot(table(vec6))

vec7 = vec^7
plot(table(vec7))

vec8 = vec^8
plot(table(vec8))

vec9 = vec^9
plot(table(vec9))

vec10 = vec^10
plot(table(vec10))

vec11 = vec^11
plot(table(vec11))

vec12 = vec^12
plot(table(vec12))

vec13 = vec^13
plot(table(vec13))

vec14 = vec^14
plot(table(vec14))


# k=4
U4 = vec+vec2
plot(table(U4))

# k=5
U5 = vec+vec2+vec3
plot(table(U5))

# k=6
U6 = vec+vec2+vec3+vec4
plot(table(U6))

# k=7
U7 = vec+vec2+vec3+vec4+vec5
plot(table(U7))

# k=8
U8 = vec+vec2+vec3+vec4+vec5+vec6
plot(table(U8))

# k=9
U9 = vec+vec2+vec3+vec4+vec5+vec6+vec7
plot(table(U9))

# k=10
U10 = vec+vec2+vec3+vec4+vec5+vec6+vec7+vec8
plot(table(U10))

