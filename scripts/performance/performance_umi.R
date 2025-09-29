
### PERFORMANCE FOR MAGIC ###

set.seed(1)

DE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TPs.txt",header = FALSE)

notDE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TNs.txt",header = FALSE)

sigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/Magic/magic.csv',header = FALSE)

NsigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/Magic/Nmagic.csv',header = FALSE)


tp <- nrow(intersect(DE,sigDE))
fp <- nrow(intersect(notDE,sigDE)) 
tn <- nrow(intersect(notDE,NsigDE))
fn <- nrow(intersect(DE,NsigDE)) 


Accuracy = (tp+tn)/(tp+tn+fn+fp)
Accuracy

Precision = tp/(tp+fp)
Precision

Recall = tp/(tp+fn)
Recall

F1 = Precision*Recall
F2 = Precision + Recall
F3 = F1/F2
F1_Score = 2*F3
F1_Score

Specificity = tn/(tn+fp)
Specificity

FDR = fp/(fp+tp)
FDR



### PERFORMANCE FOR SAVER ###

set.seed(1)

DE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TPs.txt",header = FALSE)

notDE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TNs.txt",header = FALSE)

sigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/saver/saver.csv',header = FALSE)

NsigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/saver/Nsaver.csv',header = FALSE)


tp <- nrow(intersect(DE,sigDE))
fp <- nrow(intersect(notDE,sigDE)) 
tn <- nrow(intersect(notDE,NsigDE))
fn <- nrow(intersect(DE,NsigDE)) 

Accuracy = (tp+tn)/(tp+tn+fn+fp)
Accuracy

Precision = tp/(tp+fp)
Precision

Recall = tp/(tp+fn)
Recall

F1 = Precision*Recall
F2 = Precision + Recall
F3 = F1/F2
F1_Score = 2*F3
F1_Score

Specificity = tn/(tn+fp)
Specificity

FDR = fp/(fp+tp)
FDR


### PERFORMANCE FOR BIOR 2.6 ###

set.seed(1)

DE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TPs.txt",header = FALSE)

notDE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TNs.txt",header = FALSE)

sigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/bior2.6/bior26.csv',header = FALSE)

NsigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/bior2.6/Nbior26.csv',header = FALSE)


tp <- nrow(intersect(DE,sigDE))
fp <- nrow(intersect(notDE,sigDE)) 
tn <- nrow(intersect(notDE,NsigDE))
fn <- nrow(intersect(DE,NsigDE)) 


Accuracy = (tp+tn)/(tp+tn+fn+fp)
Accuracy

Precision = tp/(tp+fp)
Precision

Recall = tp/(tp+fn)
Recall

F1 = Precision*Recall
F2 = Precision + Recall
F3 = F1/F2
F1_Score = 2*F3
F1_Score

Specificity = tn/(tn+fp)
Specificity

FDR = fp/(fp+tp)
FDR


### PERFORMANCE FOR DB6 ###

set.seed(1)

DE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TPs.txt",header = FALSE)

notDE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TNs.txt",header = FALSE)

sigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/db6/db6.csv',header = FALSE)

NsigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/db6/Ndb6.csv',header = FALSE)


tp <- nrow(intersect(DE,sigDE))
fp <- nrow(intersect(notDE,sigDE)) 
tn <- nrow(intersect(notDE,NsigDE))
fn <- nrow(intersect(DE,NsigDE)) 


Accuracy = (tp+tn)/(tp+tn+fn+fp)
Accuracy

Precision = tp/(tp+fp)
Precision

Recall = tp/(tp+fn)
Recall

F1 = Precision*Recall
F2 = Precision + Recall
F3 = F1/F2
F1_Score = 2*F3
F1_Score

Specificity = tn/(tn+fp)
Specificity

FDR = fp/(fp+tp)
FDR


### PERFORMANCE FOR ENHANCE ###

set.seed(1)

DE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TPs.txt",header = FALSE)

notDE <- read.csv("/home/irenekay/Desktop/Isabel/Real/UMI/truthset/TNs.txt",header = FALSE)

sigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/enhance/enhance.tsv')
colnames(sigDE)[1]  <- "V1" 

NsigDE = read.csv('/home/irenekay/Desktop/Isabel/Real/UMI/enhance/Nenhance.tsv')
colnames(NsigDE)[1]  <- "V1"

tp <- nrow(intersect(DE,sigDE))
fp <- nrow(intersect(notDE,sigDE)) 
tn <- nrow(intersect(notDE,NsigDE))
fn <- nrow(intersect(DE,NsigDE)) 


Accuracy = (tp+tn)/(tp+tn+fn+fp)
Accuracy

Precision = tp/(tp+fp)
Precision

Recall = tp/(tp+fn)
Recall

F1 = Precision*Recall
F2 = Precision + Recall
F3 = F1/F2
F1_Score = 2*F3
F1_Score

Specificity = tn/(tn+fp)
Specificity

FDR = fp/(fp+tp)
FDR

