# Experiment data analyses
Project repository for the functions used to fit DIVERCE experiment data.

# Files
1. gause_wrapper_fixed.R: modified version of `gauseR::gause_wrapper()` which allows you to keep some parameters fixed/constant (assuming that they have already been determined) using the `offset` argument in `lm()`. Also should allow you to optimise the fit in log-space or regular-space, according to your needs.

# Workflow
1. Try stuff
2. When things are wrong, either submit an [issue](https://github.com/nonsciencemark/data-analysis-experiments/issues) or fix it and submit a pull request
3. ???
4. Profit

Use `git clone https://github.com/nonsciencemark/data-analysis-experiments` and then `git checkout -b your-branch-name` if you want to create a new branch and work on that. Then you can make a merge request or whatever it's called.

