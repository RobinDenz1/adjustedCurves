# Contributing to adjustedCurves

First of all, thanks for considering contributing to `adjustedCurves`!

`adjustedCurves` is an open source project, maintained by people who
care. We are not directly funded to do so.

## How you can contribute

There are several ways you can contribute to this project. If you want
to know more about why and how to contribute to open source projects
like this one, see this [Open Source
Guide](https://opensource.guide/how-to-contribute/).

### Share the love

Think `adjustedCurves` is useful? Let others discover it, by telling
them in person, via Twitter or a blog post.

Using `adjustedCurves` for a paper you are writing? Consider citing it.

### Ask a question

Using `adjustedCurves` and got stuck? Browse the documentation and
vignettes to see if you can find a solution. Still stuck? Post your
question as an [issue on
GitHub](https://github.com/RobinDenz1/adjustedCurves/issues/new). While
we cannot offer user support, we’ll try to do our best to address it, as
questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by
\[email\]\[[mailto:email](mailto:email)\].

### Propose an idea

Have an idea for a new `adjustedCurves` feature? Take a look at the
documentation and [issue
list](https://github.com/RobinDenz1/adjustedCurves/issues) to see if it
isn’t included or suggested yet. If not, suggest your idea as an [issue
on GitHub](https://github.com/RobinDenz1/adjustedCurves/issues/new).
While we can’t promise to implement your idea, it helps to:

- Explain in detail how it would work.
- Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug

Using `adjustedCurves` and discovered a bug? That’s annoying! Don’t let
others have the same experience and report it as an [issue on
GitHub](https://github.com/RobinDenz1/adjustedCurves/issues/new) so we
can fix it. A good bug report makes it easier for us to do so, so please
include:

- Your operating system name and version (e.g. Mac OS 10.13.6).
- Any details about your local setup that might be helpful in
  troubleshooting.
- Detailed steps to reproduce the bug.

It is very important that a minimal reproducible example is included.
Otherwise it will be very hard for us to fix the problem.

### Improve the documentation

Noticed a typo on the website? Think a function could use a better
example? Good documentation makes all the difference, so your help to
improve it is very welcome!

### Contribute code

Care to fix bugs or implement new functionality for `adjustedCurves`?
Awesome! Have a look at the [issue
list](https://github.com/RobinDenz1/adjustedCurves/issues) and leave a
comment on the things you want to work on. See also the development
guidelines below.

## Development guidelines

Before submitting your own code it is recommended communicate with the
package maintainer about your proposal. This can be done via github or
via e-mail.

In general, we try to follow the [GitHub
flow](https://guides.github.com/introduction/flow/) for development.

1.  Fork [this repo](https://github.com/RobinDenz1/adjustedCurves) and
    clone it to your computer. To learn more about this process, see
    [this guide](https://guides.github.com/activities/forking/).
2.  If you have forked and cloned the project before and it has been a
    while since you worked on it, [pull changes from the original
    repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/)
    to your clone by using `git pull upstream master`.
3.  Open the RStudio project file (`.Rproj`).
4.  Make your changes:
    - Write your code.
    - Test your code (bonus points for adding unit tests).
    - Document your code (see function documentation above).
    - Check your code with `devtools::check()` and aim for 0 errors and
      warnings.
5.  Commit and push your changes.
6.  Submit a [pull
    request](https://guides.github.com/activities/forking/#making-a-pull-request).
