context("training")

test_that("Sucessful training from a dummy model", {
 data("parenthesis")
 maxlen <- 30
 batch.size <- 10
 preprocessed <- preprocessSemiRedundant(substr(parenthesis, 1, 100),
                                                maxlen = maxlen)
 
 expect_error(trainNetwork(""))
 expect_error(trainNetwork(dataset = preprocessed, maxlen = 0))
 expect_error(trainNetwork(dataset = preprocessed, maxlen = ""))
 expect_error(trainNetwork(dataset = preprocessed, dropout.rate = ""))
 expect_error(trainNetwork(dataset = preprocessed, dropout.rate = 0))
 expect_error(trainNetwork(dataset = preprocessed, dropout.rate = 1))
 expect_error(trainNetwork(dataset = preprocessed, layer.size = ""))
 expect_error(trainNetwork(dataset = preprocessed, layer.size = 1))
 expect_error(trainNetwork(dataset = preprocessed, layer_lstm = ""))
 expect_error(trainNetwork(dataset = preprocessed, layer_lstm  = 1))
 expect_error(trainNetwork(dataset = preprocessed, batch.size = ""))
 expect_error(trainNetwork(dataset = preprocessed, batch.size  = 1))
 expect_error(trainNetwork(dataset = "", path = ""))

 trainedNetwork <- trainNetwork(dataset = preprocessed,
                                vocabulary.size = 7,
                                batch.size = batch.size,
                                maxlen = maxlen,
                                layers.lstm = 2,
                                layer.size = 10,
                                epochs = 1)

 expect_type(trainedNetwork, "list")
 expect_equal(length(trainedNetwork),4)
 expect_type(trainedNetwork[1], "list")
 expect_equal(length(trainedNetwork[[1]]),7)
 expect_type(trainedNetwork[2], "list")
 expect_equal(length(trainedNetwork[[2]]),2)

 expect_type(trainedNetwork[[1]][["batch_size"]],"integer")
 expect_equal(trainedNetwork[[1]][["batch_size"]],10)
 expect_type(trainedNetwork[[1]][["epochs"]],"integer")
 expect_equal(trainedNetwork[[1]][["epochs"]],1)
 expect_type(trainedNetwork[[1]][["metrics"]],"character")
 expect_equal(trainedNetwork[[1]][["metrics"]],c("loss","val_loss", "acc", "val_acc"))

 expect_type(trainedNetwork[[2]][["loss"]],"double")
 expect_type(trainedNetwork[[2]][["val_loss"]],"double")
})
