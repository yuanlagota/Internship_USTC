def query(prompt, dtype=None, default=None, choices=None):
    """
    Ask for user input.

    :prompt:  default message before the input
    :dtype:   type to cast user intput in (from `default` if available)
    :default: default values if no input is given
    :choices: restrict user input to these choices
    """

    # update prompt: show list of choices
    if choices:
        prompt += " ({}".format(choices[0])
        for c in choices[1:]:
            prompt += ", {}".format(c)
        prompt += ")"

    # update prompt: show default value
    if default is not None:
        prompt += " [{}]".format(default)
        if choices  and  not default in choices:
            raise(ValueError("default value is not in list of choices"))

    prompt += ": "


    # read user input
    while True:
        s = input(prompt)
        # no input -> set default value (if available, or continue otherwise)
        if not s:
            if default is None: continue
            s = default
            break

        # convert type (if dtype is not given, then set it from default value)
        if not dtype  and  default is not None:
            dtype = type(default)
        if dtype:
            try:
                s = dtype(s)
            except Exception as error:
                print(error)
                continue

        # check if input is in list of choices
        if choices  and  not s in choices:
            print("ERROR: {} is not a valid choice!".format(s))
            continue
        break
    return s
