setwd("~/Downloads/DL")

# Read in the training files

usf1.gm.dnase.train.NY <- read.csv("usf1.gm12878.dnase.train.csv",sep=",",header=T)
usf1.gm.dnase.train.Y  <- read.csv("usf1.gm12878.Y.train.csv",sep=",",header=T)

# Count the binding outcomes
table(usf1.gm.dnase.train.Y$x)

# Create a single train data frame with binding outcomes
usf1.gm.dnase.train <- cbind(usf1.gm.dnase.train.NY,usf1.gm.dnase.train.Y)
nrow(usf1.gm.dnase.train)

# Read in the test files
usf1.gm.dnase.test.NY <- read.csv("usf1.gm12878.dnase.test.csv",sep=",",header=T)
usf1.gm.dnase.test.Y  <- read.csv("usf1.gm12878.Y.test.csv",sep=",",header=T)

# Count the binding outcomes
table(usf1.gm.dnase.test.Y$x)

# Create a single test data frame with binding outcomes
usf1.gm.dnase.test <- cbind(usf1.gm.dnase.test.NY,usf1.gm.dnase.test.Y)

# Now combine both test and training into single data frame

usf1.gm.dnase <- rbind(usf1.gm.dnase.train,usf1.gm.dnase.test)
nrow(usf1.gm.dnase)

# Mix up the occurrences of binding site 0/1 to guarantee mixture of 0 and 1s in
# the test and training sets

usf1.gm.dnase.shuff <- usf1.gm.dnase[sample(1:nrow(usf1.gm.dnase)),]
usf1.gm.dnase.shuff$x <- factor(usf1.gm.dnase.shuff$x,levels=0:1,labels=c("FALSE","TRUE"))

# Okay load up the caret package - we could just use glm straight away
# but let's see what caret has

# http://www.utstat.toronto.edu/~brunner/oldclass/appliedf11/handouts/2101f11StepwiseLogisticR.pdf
library(caret)

Train <- createDataPartition(usf1.gm.dnase.shuff$x, p=0.6, list=FALSE)
training <- usf1.gm.dnase.shuff[  Train, ]
testing  <- usf1.gm.dnase.shuff[ -Train, ]

# Let's work with the caret package

full_fit_caret <- train(x ~ .,  data=training, method="glm", family="binomial")
summary(full_fit_caret)
plot(varImp(full_fit_caret))

preds <- predict(full_fit_caret, testing)
probs <- ifelse(preds >= 0.5,TRUE,FALSE)
confusionMatrix(testing$x,probs)


full_fit <- glm(x ~ ., data=training, family="binomial")
nothing  <- glm(x ~ 1, data=training, family="binomial")

# Let's look at the full model first

full_fit_probs <- predict(full_fit,testing,type="response")
full_fit_preds <- ifelse(full_fit_probs > 0.5,TRUE,FALSE)
myt <- table(full_fit_preds,testing$x)
mym <- mean(full_fit_preds==testing$x)
myt
mym

# But most of the coefficients are not significant


full_fit2 <- glm(x ~ V2 + V10 + V11 + V13 + V20, data=training,family="binomial")

full_fit2_probs <- predict(full_fit2,testing,type="response")
full_fit2_preds <- ifelse(full_fit2_probs > 0.5,TRUE,FALSE)
myt <- table(full_fit_preds,testing$x)
mym <- mean(full_fit_preds==testing$x)
myt
mym

library(ROCR)
pred <- prediction(full_fit2_probs,testing$x)
perf <- performance(pred,"tpr","fpr")
plot(perf)

# Stepwise is possible here but let's hold off for now
backwards <- step(full_fit)

forwards  <- step(nothing,scope=list(lower=formula(nothing),upper=formula(full_fit)),direction="forward")

library(pROC)

mythres <- function(mod, data, thres=0.5) {
  myt <- table(predict(mod,data) > thres,data$x)
  return(myt)
}

mythres(backwards,testing,0.5)

roc(testing$x,predict(forwards,testing)) -> rroc
plot(rroc,type="S",print.thres=.9)

# rpart

cfit <- randomForest(x~V11 + V10 + V17 + V12 + V2 + V9 + V7 + V20 + V19,data=training,
                     + importance=TRUE,ntree=2000)

cfit <- randomForest(x~.,data=training,
                     importance=TRUE,ntree=2000)

myglm <- glm(x ~ V14 + V8 + V7 + V14 + V13 + V6 + V12 + V9 + V11 + V10,
             data=training,family="binomial")

myglm.probs <- predict(myglm,testing,type="response")
myglm.preds <- ifelse(myglm.probs > 0.5,TRUE,FALSE)
myt <- table(myglm.preds,testing$x)
mym <- mean(myglm.preds==testing$x)
