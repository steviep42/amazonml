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

library(caret)

Train <- createDataPartition(usf1.gm.dnase.shuff$x, p=0.6, list=FALSE)
training <- usf1.gm.dnase.shuff[  Train, ]
testing  <- usf1.gm.dnase.shuff[ -Train, ]

myglm <- glm(x ~ .,data=training,family="binomial")
summary(myglm)

df <- varImp(myglm)

# Sort by Overall variable importance

df <- df[order(df$Overall),,drop=FALSE]
dotchart(df$Overall,labels=rownames(df))
dotchart(df$Overall,labels=rownames(df),pch=19,col="blue",main="Variable Importance Plot")

# Let's create some predictions and then generate a confusion matrix
# We write a little helper function to let us vary the threshold used
# to compute the predictions of TRUE or FALSE (1 or 0)

myglm.probs <- predict(myglm, testing, type="response")
myglm.preds <- function(t) ifelse(myglm.probs > t, TRUE, FALSE)

confusionMatrix(myglm.preds(0.5),testing$x)

#Let's just get the accuracy 

confusionMatrix(myglm.preds(0.5),testing$x)$overall[1]

# We could run this for a variety of thresholds

acc <- sapply(seq(0.2,0.8,.05),
              function(t) confusionMatrix(myglm.preds(t),testing$x)$overall[1] )

# Name each element value according to the respective threshold 
names(acc) <- seq(0.2,0.8,.05)

(round(acc,3))

# Which treshold corresponds to the highest degree of accuracy ? 
which(acc==max(acc))

# Let's draw the ROC curve

library(ROCR)

preds <- prediction(myglm.probs, testing$x)
perf  <- performance(preds, measure = "tpr", x.measure = "fpr") 
plot(perf, colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7))

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
title(paste("ROC Curve with AUC =",round(auc,3)))


# Cross validation

mod_fit <- train(x ~ .,data=training,method="glm",
                 family="binomial",trControl=ctrl,
                 tuneLength=5)

fitted.vals <- mod_fit$finalModel$fitted.values
fitpredt <- function(t) ifelse(fitted.vals > t , TRUE, FALSE)
confusionMatrix(fitpredt(0.5), training$x)

summary(mod_fit)

# Random forests

library(randomForest)

forest_fit <- randomForest(x~.,data=training, importance=TRUE,ntree=2000)

varImpPlot(cfit,pch=19,col="blue")

predictions <- predict(forest_fit, testing)

confusionMatrix(predictions,testing$x)

