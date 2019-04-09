# Awesome
## linux editing tricks
### linux terminal sed and tr example

Pasted from linux terminal:

```
-bash-4.2$ echo "Autism_123, Typical_456"
Autism_123, Typical_456
-bash-4.2$ echo "Autism_123, Typical_456"|sed 's/,//g'
Autism_123 Typical_456
-bash-4.2$ echo "Autism_123, Typical_456"|sed 's/,//g'|tr ' ' '\n'
Autism_123
Typical_456
-bash-4.2$ echo "Autism_123, Typical_456"|sed 's/,//g'|tr ' ' '\n'|sed 's#\(.*\)#autism_\1_HD <- estimateDivergence(Ref_13#'
autism_Autism_123_HD <- estimateDivergence(Ref_13
autism_Typical_456_HD <- estimateDivergence(Ref_13
```
