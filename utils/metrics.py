import tensorflow as tf
from sklearn.metrics import matthews_corrcoef, f1_score, roc_auc_score
from tensorflow.keras.losses import binary_crossentropy


def masked_softmax_cross_entropy(preds, labels, mask):
    """
    Softmax cross-entropy loss with masking.
    """
    # loss = tf.nn.softmax_cross_entropy_with_logits(logits=preds, labels=labels)
    # mask = tf.cast(mask, dtype=tf.float32)
    # mask /= tf.reduce_mean(mask)
    # loss *= mask
    # return tf.reduce_mean(loss)
    curr_labels = tf.boolean_mask(labels, mask)
    curr_preds = tf.nn.sigmoid(tf.boolean_mask(preds, mask))
    loss2 = tf.compat.v1.losses.log_loss(curr_labels, curr_preds)
    return loss2



def masked_accuracy(preds, labels, mask):
    """
    Accuracy with masking.
    """
    correct_prediction = tf.equal(tf.argmax(preds, 1), tf.argmax(labels, 1))
    accuracy_all = tf.cast(correct_prediction, tf.float32)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    accuracy_all *= mask
    return tf.reduce_mean(accuracy_all)


def masked_f1_score(preds, labels, mask):
    """
    F1 score with masking.
    """
    curr_labels = tf.boolean_mask(labels, mask)
    curr_preds = tf.nn.softmax(tf.boolean_mask(preds, mask), axis=1)

    curr_preds = tf.math.greater_equal(curr_preds, tf.constant(0.5, dtype='float32'))
    return f1_score(y_true=curr_labels.numpy(), y_pred=curr_preds.numpy())


def masked_matthews_corrcoef(preds, labels, mask):
    """
    MCC score with masking.
    """
    # curr_labels = tf.boolean_mask(labels, mask)[:, 0]
    # curr_preds = tf.nn.softmax(tf.boolean_mask(preds, mask), axis=1)[:, 0]
    #
    # curr_preds = tf.math.greater_equal(curr_preds, tf.constant(0.5, dtype='float32'))
    # curr_preds = tf.cast(curr_preds, dtype=tf.float32)

    curr_labels = tf.boolean_mask(labels, mask)
    curr_preds = tf.boolean_mask(preds, mask)
    curr_preds = tf.nn.sigmoid(curr_preds)
    curr_preds = tf.math.greater_equal(curr_preds, tf.constant(0.5, dtype='float32'))
    curr_preds = tf.cast(curr_preds, dtype=tf.float32)
    return matthews_corrcoef(y_true=curr_labels.numpy(), y_pred=curr_preds.numpy())


def masked_roc_auc(preds, labels, mask):
    """
    AUROC score with masking.
    """
    curr_labels = tf.boolean_mask(labels, mask)
    curr_preds = tf.boolean_mask(preds, mask)
    # curr_preds = tf.nn.sigmoid(curr_preds)
    return roc_auc_score(y_true=curr_labels.numpy().flatten(), y_score=curr_preds.numpy().flatten())


def masked_binary_cross_entropy(preds, labels, mask):
    """
    binary cross-entropy loss with masking.
    """
    loss = binary_crossentropy(y_true=labels, y_pred=preds, from_logits=True)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    loss *= mask
    return tf.reduce_mean(loss)

