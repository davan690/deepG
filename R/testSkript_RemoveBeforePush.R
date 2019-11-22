library(keras)
hist <- trainNetwork(path = "/home/rmreches/resumeTraining/parenthesis_as_fasta/train",
             path.val = "/home/rmreches/resumeTraining/parenthesis_as_fasta/validation",
             checkpoint_path = "/home/rmreches/resumeTraining/checkpoints",
             run.name = "model_3_layers",
             maxlen = 50,
             dropout = 0.2,
             recurrent_dropout = 0.2,
             layer.size = 128,
             batch.size = 32,
             layers.lstm = 2,
             vocabulary.size = 7,
             epochs = 6,
             max.queue.size = 10,
             steps.per.epoch = 5,
             step = 25,
             seqStart = "",
             seqEnd= "",
             withinFile = "",
             vocabulary = c("a","c","g", "t", "m", "n","r"),
             tensorboard.log = "/home/rmreches/resumeTraining/tb")

tensorboard("/home/rmreches/resumeTraining/tb")

resumeTraining(model_path = "/home/rmreches/resumeTraining/checkpoints/runPar1_checkpoints/Ep.010-1.47.hdf5",
               path = "/home/rmreches/resumeTraining/parenthesis_as_fasta/train",
               path.val = "/home/rmreches/resumeTraining/parenthesis_as_fasta/validation",
               checkpoint_path = "/home/rmreches/resumeTraining/checkpoints",
               run.name = "runPar1_cont_from_Ep5",
               batch.size = 32,
               epochs = 15,
               seqStart = "",
               seqEnd= "",
               withinFile = "",
               max.queue.size = 10,
               steps.per.epoch = 50,
               step = 25,
               #initial_epoch = 14,
               vocabulary = c("a", "c", "g", "t", "m", "n", "r"),
               tensorboard.log = "/home/rmreches/resumeTraining/tb",
               compile = TRUE)

model <- keras::load_model_hdf5(model_path, compile = TRUE)

gen <- fastaFileGenerator(corpus.dir = "/home/rmreches/resumeTraining/parenthesis_as_fasta/train",
                          format = "fasta",
                          batch.size = 1,
                          maxlen = 50,
                          max_iter = 20,
                          seqStart = "",
                          seqEnd= "",
                          withinFile = "",
                          vocabulary = c("a", "c", "g", "t", "m", "n", "r"),
                          verbose = FALSE,
                          randomFiles = FALSE,
                          step = 1, 
                          showWarnings = FALSE)

a <- gen()
d <- a[[1]] + 0.0
which.max(predict(model, d))


